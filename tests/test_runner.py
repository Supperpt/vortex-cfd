"""Tests for the checkMesh quality gate (pure string parsing, no OpenFOAM)."""
import pytest

import vortex_cfd.runner as runner
from vortex_cfd.runner import _check_mesh_quality


def _checkmesh_output(nonortho: float, skewness: float) -> str:
    """Reproduce the relevant checkMesh lines in their real ESI v2406 format."""
    return (
        "Checking geometry...\n"
        f"    Mesh non-orthogonality Max: {nonortho} average: 12.28\n"
        "    Non-orthogonality check OK.\n"
        "    Face pyramids OK.\n"
        f" ***Max skewness = {skewness}, 29 highly skew faces detected\n"
        "    Coupled point location match (average 0) OK.\n"
    )


def test_good_mesh_passes(capsys):
    # AA_011: skewness 6.66, non-ortho 69.96 — both under the gates.
    _check_mesh_quality(_checkmesh_output(69.96, 6.66))
    assert "PASSED" in capsys.readouterr().out


def test_skewness_under_20_passes():
    # Anything up to 20 (OpenFOAM maxBoundarySkewness) is acceptable.
    _check_mesh_quality(_checkmesh_output(50.0, 19.9))  # no SystemExit


def test_skewness_over_20_aborts():
    with pytest.raises(SystemExit):
        _check_mesh_quality(_checkmesh_output(50.0, 25.0))


def test_nonortho_over_70_aborts():
    # Regression: the non-ortho regex must match checkMesh's real output format.
    with pytest.raises(SystemExit):
        _check_mesh_quality(_checkmesh_output(75.0, 3.0))


def test_old_skewness_threshold_no_longer_aborts():
    # 6.66 used to trip the old maxSkewness>4 gate (BUG-009); must pass now.
    _check_mesh_quality(_checkmesh_output(69.96, 6.66))  # no SystemExit


# --- run_log.md (run-step logging) ------------------------------------------

def test_log_append_is_noop_without_init(tmp_path, monkeypatch):
    # When no log is initialised, _log_append must do nothing (and not raise).
    monkeypatch.setattr(runner, "_RUN_LOG", None)
    runner._log_append("## should not be written\n")  # no error, nothing created
    assert not list(tmp_path.iterdir())


def test_log_init_and_append_write_markdown(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "_RUN_LOG", None)
    runner._log_init(tmp_path, "full pipeline")
    runner._log_append("\n## blockMesh — OK\n\n```\nmesh stats...\n```\n")
    text = (tmp_path / "run_log.md").read_text()
    assert "# vortex-cfd run log" in text
    assert "full pipeline" in text
    assert "## blockMesh — OK" in text
