"""Tests for vortex_cfd.case_builder geometry helpers and full case generation."""

import json

import numpy as np
import pyvista as pv
import pytest

from vortex_cfd.case_builder import (
    _background_cell_counts,
    _bbox_with_buffer,
    _inlet_area,
    _location_in_mesh,
    _waveform_table,
    build_case,
)
from vortex_cfd.waveform import T_CYCLE, load_waveform


# ---------------------------------------------------------------------------
# _bbox_with_buffer
# ---------------------------------------------------------------------------

class TestBboxWithBuffer:
    def test_buffered_bbox_is_larger(self, scaled_stls_m):
        wall = pv.read(str(scaled_stls_m["wall"]))
        b = wall.bounds
        bbox = _bbox_with_buffer(scaled_stls_m["wall"], buffer=0.20)
        assert bbox["xmin"] < b[0]
        assert bbox["xmax"] > b[1]
        assert bbox["ymin"] < b[2]
        assert bbox["ymax"] > b[3]
        assert bbox["zmin"] < b[4]
        assert bbox["zmax"] > b[5]

    def test_buffer_size_is_20_percent(self, scaled_stls_m):
        wall = pv.read(str(scaled_stls_m["wall"]))
        b = wall.bounds
        dx_orig = b[1] - b[0]
        bbox = _bbox_with_buffer(scaled_stls_m["wall"], buffer=0.20)
        dx_new = bbox["xmax"] - bbox["xmin"]
        assert dx_new == pytest.approx(dx_orig * 1.40, rel=1e-4)

    def test_returns_all_six_keys(self, scaled_stls_m):
        bbox = _bbox_with_buffer(scaled_stls_m["wall"])
        assert set(bbox.keys()) == {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"}

    def test_bbox_is_valid(self, scaled_stls_m):
        bbox = _bbox_with_buffer(scaled_stls_m["wall"])
        assert bbox["xmin"] < bbox["xmax"]
        assert bbox["ymin"] < bbox["ymax"]
        assert bbox["zmin"] < bbox["zmax"]


# ---------------------------------------------------------------------------
# _background_cell_counts
# ---------------------------------------------------------------------------

class TestBackgroundCellCounts:
    def test_minimum_cells_per_direction(self):
        # Very small geometry → floor at 10 cells
        tiny = {"xmin": 0, "xmax": 0.001, "ymin": 0, "ymax": 0.001, "zmin": 0, "zmax": 0.001}
        counts = _background_cell_counts(tiny)
        for v in counts.values():
            assert v >= 10

    def test_maximum_cells_per_direction(self):
        # Very large geometry → cap at 60 cells
        huge = {"xmin": 0, "xmax": 1.0, "ymin": 0, "ymax": 1.0, "zmin": 0, "zmax": 1.0}
        counts = _background_cell_counts(huge)
        for v in counts.values():
            assert v <= 60

    def test_typical_geometry_reasonable_counts(self):
        # 30 mm vessel field → expect ~15 cells with 2 mm target
        typical = {"xmin": -0.012, "xmax": 0.018,
                   "ymin": -0.010, "ymax": 0.020,
                   "zmin": -0.005, "zmax": 0.025}
        counts = _background_cell_counts(typical)
        for v in counts.values():
            assert 10 <= v <= 60

    def test_returns_nx_ny_nz(self):
        bbox = {"xmin": 0, "xmax": 0.03, "ymin": 0, "ymax": 0.02, "zmin": 0, "zmax": 0.04}
        counts = _background_cell_counts(bbox)
        assert set(counts.keys()) == {"nx", "ny", "nz"}


# ---------------------------------------------------------------------------
# _inlet_area
# ---------------------------------------------------------------------------

class TestInletArea:
    def test_positive(self, scaled_stls_m):
        area = _inlet_area(scaled_stls_m["inlet"])
        assert area > 0

    def test_disc_area_close_to_pi_r_squared(self, scaled_stls_m):
        # Inlet disc was created with radius 0.003 m → π×r² ≈ 2.827e-5 m²
        area = _inlet_area(scaled_stls_m["inlet"])
        expected = np.pi * 0.003 ** 2
        # 5 % tolerance — triangulated disc underestimates slightly
        assert area == pytest.approx(expected, rel=0.05)


# ---------------------------------------------------------------------------
# _location_in_mesh
# ---------------------------------------------------------------------------

class TestLocationInMesh:
    def test_returns_three_floats(self, scaled_stls_m):
        loc = _location_in_mesh(scaled_stls_m["wall"])
        assert len(loc) == 3
        assert all(isinstance(v, float) for v in loc)

    def test_location_inside_wall_bbox(self, scaled_stls_m):
        wall = pv.read(str(scaled_stls_m["wall"]))
        b = wall.bounds
        loc = _location_in_mesh(scaled_stls_m["wall"])
        assert b[0] <= loc[0] <= b[1]
        assert b[2] <= loc[1] <= b[3]
        assert b[4] <= loc[2] <= b[5]


# ---------------------------------------------------------------------------
# _waveform_table
# ---------------------------------------------------------------------------

class TestWaveformTable:
    def test_length_matches_waveform(self):
        wf = load_waveform(None)
        table = _waveform_table(wf, mean_velocity=0.4, inlet_area=1e-4)
        assert len(table) == len(wf)

    def test_mean_flow_rate_equals_u_times_a(self):
        wf = load_waveform(None)
        U, A = 0.4, 1.5e-5
        table = _waveform_table(wf, mean_velocity=U, inlet_area=A)
        Q_vals = [q for _, q in table]
        assert np.mean(Q_vals) == pytest.approx(U * A, rel=1e-4)

    def test_times_are_in_seconds(self):
        wf = load_waveform(None)
        table = _waveform_table(wf, mean_velocity=0.4, inlet_area=1e-4)
        t_last = table[-1][0]
        # Last time should be close to (but less than) one full cycle
        assert 0 < t_last < T_CYCLE

    def test_all_flow_rates_positive(self):
        wf = load_waveform(None)
        table = _waveform_table(wf, mean_velocity=0.4, inlet_area=1e-4)
        assert all(q > 0 for _, q in table)

    def test_flow_rate_scales_linearly_with_velocity(self):
        wf = load_waveform(None)
        A = 1e-4
        t1 = _waveform_table(wf, mean_velocity=0.3, inlet_area=A)
        t2 = _waveform_table(wf, mean_velocity=0.6, inlet_area=A)
        for (_, q1), (_, q2) in zip(t1, t2):
            assert q2 == pytest.approx(q1 * 2.0, rel=1e-6)


# ---------------------------------------------------------------------------
# build_case — full case generation
# ---------------------------------------------------------------------------

REQUIRED_FILES = [
    "0/U",
    "0/p",
    "constant/transportProperties",
    "constant/turbulenceProperties",
    "system/controlDict",
    "system/fvSchemes",
    "system/fvSolution",
    "system/snappyHexMeshDict",
    "system/decomposeParDict",
    "system/meshQualityDict",
    "system/surfaceFeatureExtractDict",
    "system/blockMeshDict",
    "Allrun",
    "patch_labels.json",
]


@pytest.fixture
def built_case(scaled_stls_m, tmp_path):
    wf = load_waveform(None)
    return build_case(
        scaled_stls=scaled_stls_m,
        labels={},
        cycles=1,
        mean_velocity=0.4,
        waveform=wf,
        cores=2,
        out_dir=str(tmp_path),
    )


class TestBuildCase:
    def test_case_directory_created(self, built_case):
        assert built_case.exists()
        assert built_case.is_dir()

    def test_case_name_starts_with_case(self, built_case):
        assert built_case.name.startswith("case_")

    @pytest.mark.parametrize("rel_path", REQUIRED_FILES)
    def test_required_file_exists(self, built_case, rel_path):
        assert (built_case / rel_path).exists(), f"Missing: {rel_path}"

    @pytest.mark.parametrize("rel_path", REQUIRED_FILES)
    def test_required_file_non_empty(self, built_case, rel_path):
        assert (built_case / rel_path).stat().st_size > 0, f"Empty: {rel_path}"

    def test_foam_placeholder_exists(self, built_case):
        foam_files = list(built_case.glob("*.foam"))
        assert len(foam_files) == 1

    def test_stls_copied_to_triSurface(self, built_case):
        ts = built_case / "constant" / "triSurface"
        assert (ts / "wall.stl").exists()
        assert (ts / "inlet.stl").exists()
        assert (ts / "outlet_0.stl").exists()

    def test_patch_labels_json_valid(self, built_case):
        meta = json.loads((built_case / "patch_labels.json").read_text())
        assert isinstance(meta, dict)

    def test_patch_labels_contains_wall(self, built_case):
        meta = json.loads((built_case / "patch_labels.json").read_text())
        assert "wall" in meta.values()

    def test_patch_labels_contains_inlet(self, built_case):
        meta = json.loads((built_case / "patch_labels.json").read_text())
        assert "inlet" in meta.values()

    def test_patch_labels_contains_outlet(self, built_case):
        meta = json.loads((built_case / "patch_labels.json").read_text())
        assert "outlet" in meta.values()

    def test_control_dict_end_time(self, built_case):
        # 1 cycle → endTime ≈ T_CYCLE (0.857 s)
        text = (built_case / "system" / "controlDict").read_text()
        assert f"{T_CYCLE:.4f}" in text

    def test_decompose_par_cores(self, built_case):
        text = (built_case / "system" / "decomposeParDict").read_text()
        assert "2" in text  # cores=2 passed to fixture

    def test_U_contains_waveform_table(self, built_case):
        text = (built_case / "0" / "U").read_text()
        assert "flowRateInletVelocity" in text
        assert "table" in text
        assert "outOfBounds" in text

    def test_p_contains_fixed_value_at_outlet(self, built_case):
        text = (built_case / "0" / "p").read_text()
        assert "fixedValue" in text

    def test_snappy_contains_location_in_mesh(self, built_case):
        text = (built_case / "system" / "snappyHexMeshDict").read_text()
        assert "locationInMesh" in text

    def test_block_mesh_dict_has_vertices(self, built_case):
        text = (built_case / "system" / "blockMeshDict").read_text()
        assert "vertices" in text
        assert "blocks" in text
