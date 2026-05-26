"""Tests for vortex_cfd.scaling."""

import pyvista as pv
import pytest

from vortex_cfd.scaling import MM_TO_M, _any_in_mm, scale_stls


class TestDetection:
    def test_metre_geometry_not_detected_as_mm(self, stl_paths_m):
        assert _any_in_mm(stl_paths_m) is False

    def test_mm_geometry_detected(self, stl_paths_mm):
        assert _any_in_mm(stl_paths_mm) is True

    def test_single_large_file_triggers_detection(self, stl_paths_mm):
        # Even if only one file in the list is in mm, detection fires
        assert _any_in_mm([stl_paths_mm[0]]) is True


class TestScaling:
    def test_metre_geometry_not_scaled(self, stl_paths_m, labels_m):
        result = scale_stls(stl_paths_m, labels_m)
        wall = pv.read(str(result["wall"]))
        # Bounding-box max should still be < 1.0 (metres)
        assert max(abs(v) for v in wall.bounds) < 1.0

    def test_mm_geometry_scaled_to_metres(self, stl_paths_mm, labels_mm):
        result = scale_stls(stl_paths_mm, labels_mm)
        wall = pv.read(str(result["wall"]))
        # After scaling, bounding-box max must be < 1.0 (metres)
        assert max(abs(v) for v in wall.bounds) < 1.0

    def test_scaling_factor_is_correct(self, stl_paths_mm, labels_mm):
        # Wall sphere was r=8.0 mm → should become r≈0.008 m after scaling
        result = scale_stls(stl_paths_mm, labels_mm)
        wall = pv.read(str(result["wall"]))
        max_coord = max(abs(v) for v in wall.bounds)
        assert max_coord == pytest.approx(0.008, rel=0.02)

    def test_metre_geometry_coords_unchanged(self, stl_paths_m, labels_m):
        # Read original wall
        wall_orig_path = next(p for p in stl_paths_m if "wall" in p.name)
        orig_max = max(abs(v) for v in pv.read(str(wall_orig_path)).bounds)

        result = scale_stls(stl_paths_m, labels_m)
        scaled_max = max(abs(v) for v in pv.read(str(result["wall"])).bounds)
        assert scaled_max == pytest.approx(orig_max, rel=1e-4)


class TestCanonicalNames:
    def test_wall_key_present(self, stl_paths_m, labels_m):
        result = scale_stls(stl_paths_m, labels_m)
        assert "wall" in result

    def test_inlet_key_present(self, stl_paths_m, labels_m):
        result = scale_stls(stl_paths_m, labels_m)
        assert "inlet" in result

    def test_outlet_numbered_from_zero(self, stl_paths_m, labels_m):
        result = scale_stls(stl_paths_m, labels_m)
        assert "outlet_0" in result

    def test_no_original_names_in_keys(self, stl_paths_m, labels_m):
        result = scale_stls(stl_paths_m, labels_m)
        for key in result:
            assert key in {"wall", "inlet", "outlet_0", "outlet_1", "outlet_2"}

    def test_multiple_outlets_numbered_sequentially(self, tmp_path):
        """Three outlet caps should become outlet_0, outlet_1, outlet_2."""
        from tests.conftest import _sphere_stl, _disc_stl
        _sphere_stl(tmp_path / "wall.stl", radius=0.008)
        _disc_stl(tmp_path / "inlet.stl", radius=0.003)
        _disc_stl(tmp_path / "out_a.stl", radius=0.002)
        _disc_stl(tmp_path / "out_b.stl", radius=0.002)
        _disc_stl(tmp_path / "out_c.stl", radius=0.002)

        paths = sorted(tmp_path.glob("*.stl"))
        labels = {}
        for p in paths:
            if "wall" in p.name:
                labels[p] = "wall"
            elif "inlet" in p.name:
                labels[p] = "inlet"
            else:
                labels[p] = "outlet"

        result = scale_stls(paths, labels)
        assert "outlet_0" in result
        assert "outlet_1" in result
        assert "outlet_2" in result

    def test_output_files_exist(self, stl_paths_m, labels_m):
        result = scale_stls(stl_paths_m, labels_m)
        for _, path in result.items():
            assert path.exists(), f"Expected output file missing: {path}"

    def test_output_files_are_readable_stls(self, stl_paths_m, labels_m):
        result = scale_stls(stl_paths_m, labels_m)
        for canonical, path in result.items():
            mesh = pv.read(str(path))
            assert mesh.n_points > 0, f"{canonical}: empty mesh after scaling"
