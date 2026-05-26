"""Shared pytest fixtures — synthetic STL geometry, no real patient data needed."""

import numpy as np
import pytest
import pyvista as pv


# ---------------------------------------------------------------------------
# STL geometry helpers
# ---------------------------------------------------------------------------

def _sphere_stl(path, radius):
    mesh = pv.Sphere(radius=radius, theta_resolution=20, phi_resolution=20)
    mesh.save(str(path), binary=False)
    return path


def _disc_stl(path, radius, center=(0.0, 0.0, 0.0)):
    """Flat triangulated disc centred at `center` lying in the XY plane."""
    n = 32
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
    rim = np.column_stack([
        center[0] + radius * np.cos(theta),
        center[1] + radius * np.sin(theta),
        np.full(n, center[2]),
    ])
    pts = np.vstack([np.array([[*center]]), rim])
    faces = np.array([v for i in range(n) for v in [3, 0, i + 1, (i % n) + 1 + 1 - (1 if i == n - 1 else 0)]])
    # Rebuild faces cleanly
    face_list = []
    for i in range(n):
        j = (i + 1) % n
        face_list.extend([3, 0, i + 1, j + 1])
    disc = pv.PolyData(pts, np.array(face_list))
    disc.save(str(path), binary=False)
    return path


# ---------------------------------------------------------------------------
# Metre-scale fixture set (no scaling should be applied)
# ---------------------------------------------------------------------------

@pytest.fixture
def stl_dir_m(tmp_path):
    """
    Temporary directory with wall + inlet + outlet_0 STLs in metres.
    Wall: sphere r=0.008 m.  Caps: discs r=0.003 m and r=0.0025 m.
    All bounding-box values < 1.0  → scaling must NOT be triggered.
    """
    wall_path = tmp_path / "wall_surface.stl"
    inlet_path = tmp_path / "cap_00.stl"
    outlet_path = tmp_path / "cap_01.stl"
    _sphere_stl(wall_path, radius=0.008)
    _disc_stl(inlet_path, radius=0.003, center=(0.0, 0.0, 0.0))
    _disc_stl(outlet_path, radius=0.0025, center=(0.005, 0.0, 0.0))
    return tmp_path


@pytest.fixture
def stl_paths_m(stl_dir_m):
    return sorted(stl_dir_m.glob("*.stl"))


@pytest.fixture
def labels_m(stl_paths_m):
    """Map each path to a semantic label for the metre-scale set."""
    d = {}
    for p in stl_paths_m:
        if "wall" in p.name:
            d[p] = "wall"
        elif "00" in p.name:
            d[p] = "inlet"
        else:
            d[p] = "outlet"
    return d


@pytest.fixture
def scaled_stls_m(stl_paths_m, labels_m):
    """Output of scale_stls for the metre-scale set (canonical names, no scaling)."""
    from vortex_cfd.scaling import scale_stls
    return scale_stls(stl_paths_m, labels_m)


# ---------------------------------------------------------------------------
# Millimetre-scale fixture set (scaling MUST be applied)
# ---------------------------------------------------------------------------

@pytest.fixture
def stl_dir_mm(tmp_path):
    """
    Same geometry as stl_dir_m but with coordinates in millimetres.
    Wall: sphere r=8 mm.  Caps: discs r=3 mm and r=2.5 mm.
    Bounding-box max ≈ 8 > 1.0  → scaling MUST be triggered.
    """
    d = tmp_path / "mm_stls"
    d.mkdir()
    _sphere_stl(d / "wall_surface.stl", radius=8.0)
    _disc_stl(d / "cap_00.stl", radius=3.0)
    _disc_stl(d / "cap_01.stl", radius=2.5)
    return d


@pytest.fixture
def stl_paths_mm(stl_dir_mm):
    return sorted(stl_dir_mm.glob("*.stl"))


@pytest.fixture
def labels_mm(stl_paths_mm):
    d = {}
    for p in stl_paths_mm:
        if "wall" in p.name:
            d[p] = "wall"
        elif "00" in p.name:
            d[p] = "inlet"
        else:
            d[p] = "outlet"
    return d
