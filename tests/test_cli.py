"""CLI wiring tests — no OpenFOAM, no meshing. The heavy stages are patched so
we only assert how `main` forwards options (notably --cycles) to the runner.
"""
from unittest import mock

from click.testing import CliRunner

import vortex_cfd.cli as cli


_FAKE_ENV = {"version": "2406", "root": "/opt/openfoam2406", "env": {}, "bashrc": None}


def test_postprocess_only_forwards_none_cycles_when_omitted():
    """
    Regression: --postprocess-only without --cycles must pass cycles=None so
    compute_metrics derives the last cycle from the case data. Defaulting to 3
    here silently analysed the wrong window for any case not solved with exactly
    3 cycles.
    """
    with mock.patch.object(cli, "check_openfoam", return_value=_FAKE_ENV), \
         mock.patch.object(cli, "run_postprocess_only") as run_pp:
        result = CliRunner().invoke(cli.main, ["--postprocess-only", "."])

    assert result.exit_code == 0, result.output
    run_pp.assert_called_once()
    assert run_pp.call_args.kwargs["cycles"] is None


def test_postprocess_only_forwards_explicit_cycles():
    """An explicit --cycles is honoured on the post-process-only path."""
    with mock.patch.object(cli, "check_openfoam", return_value=_FAKE_ENV), \
         mock.patch.object(cli, "run_postprocess_only") as run_pp:
        result = CliRunner().invoke(cli.main, ["--postprocess-only", ".", "--cycles", "5"])

    assert result.exit_code == 0, result.output
    assert run_pp.call_args.kwargs["cycles"] == 5


def test_full_run_defaults_cycles_to_three(tmp_path):
    """A full run with --cycles omitted still solves the documented 3 cycles."""
    (tmp_path / "wall.stl").write_text("")  # only needs to exist for the glob

    with mock.patch.object(cli, "check_openfoam", return_value=_FAKE_ENV), \
         mock.patch.object(cli, "label_patches", return_value={}), \
         mock.patch.object(cli, "scale_stls", return_value={}), \
         mock.patch.object(cli, "load_waveform", return_value=None), \
         mock.patch.object(cli, "build_case", return_value=tmp_path / "case") as build, \
         mock.patch.object(cli, "run_pipeline") as run_pipe:
        result = CliRunner().invoke(
            cli.main,
            ["--stl-dir", str(tmp_path), "--mean-velocity", "0.4"],
        )

    assert result.exit_code == 0, result.output
    assert build.call_args.kwargs["cycles"] == 3
    assert run_pipe.call_args.kwargs["cycles"] == 3
