"""Tests for the checkMesh quality gate (pure string parsing, no OpenFOAM)."""
import pytest

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
