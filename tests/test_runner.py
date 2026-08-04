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


# ---------------------------------------------------------------------------
# Womersley boundaryData step
# ---------------------------------------------------------------------------

class TestWomersleyBoundaryData:
    """
    The boundaryData step must fail loudly and early. Unlike a post-processing
    diagnostic, a missing constant/boundaryData leaves the solver with an inlet
    BC it cannot read -- cheaper to abort now than hours into the solve.
    """

    def _params(self):
        import numpy as np
        from vortex_cfd.waveform import load_waveform
        return {
            "centroid": np.zeros(3),
            "normal": np.array([0.0, 0.0, 1.0]),
            "radius": 0.003,
            "mean_velocity": 0.4,
            "waveform": load_waveform(None),
            "nu": 3.3e-6,
            "cycles": 2,
        }

    def test_missing_inlet_params_aborts(self, tmp_path, capsys):
        with pytest.raises(SystemExit) as exc:
            runner._write_womersley_boundary_data(tmp_path, None)
        assert exc.value.code == 1
        assert "requires inlet geometry" in capsys.readouterr().err

    def test_unreadable_mesh_aborts(self, tmp_path, capsys):
        """No .foam file — must exit rather than propagate a raw traceback."""
        with pytest.raises(SystemExit) as exc:
            runner._write_womersley_boundary_data(tmp_path, self._params())
        assert exc.value.code == 1
        assert "could not write Womersley boundaryData" in capsys.readouterr().err

    def test_writes_boundary_data_for_the_whole_run(self, tmp_path, monkeypatch):
        import numpy as np
        from vortex_cfd.waveform import T_CYCLE

        faces = np.column_stack([
            np.linspace(0.0, 0.0025, 12), np.zeros(12), np.zeros(12)])
        monkeypatch.setattr(runner, "_inlet_face_centers", lambda *a, **k: faces)

        runner._write_womersley_boundary_data(tmp_path, self._params())

        bd = tmp_path / "constant" / "boundaryData" / "inlet"
        assert (bd / "points").exists()
        written = sorted(float(p.parent.name) for p in bd.glob("*/U"))
        # cycles=2 in _params, so the data must span two full cycles.
        assert max(written) >= 2 * T_CYCLE - 1e-9
        assert len(written) == 2 * 100 + 1
