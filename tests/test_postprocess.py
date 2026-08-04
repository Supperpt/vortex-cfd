"""
Tests for vortex_cfd.postprocess pure functions (TAWSS / OSI / summary stats).

These use synthetic numpy WSS series only — no OpenFOAM, no pyvista reader —
preserving the suite's "runs without a solver" property.  The I/O layer
(read_wss_series / compute_metrics) is exercised only by the integration smoke
test against a real solved case.
"""

import json
import textwrap

import numpy as np
import pytest

from vortex_cfd.postprocess import (
    RHO,
    _detect_patches,
    _read_surface_field_value,
    _time_weights,
    osi,
    summary_stats,
    tawss,
)


# ---------------------------------------------------------------------------
# _time_weights
# ---------------------------------------------------------------------------

class TestTimeWeights:
    def test_single_sample(self):
        assert _time_weights([0.5]).tolist() == [1.0]

    def test_uniform_spacing_sums_to_interval(self):
        t = np.linspace(0.0, 1.0, 11)
        w = _time_weights(t)
        assert w.sum() == pytest.approx(1.0)

    def test_interior_weights_uniform(self):
        # Evenly spaced → interior weights all equal to the spacing.
        t = np.linspace(0.0, 1.0, 5)  # spacing 0.25
        w = _time_weights(t)
        assert w[1] == pytest.approx(0.25)
        assert w[2] == pytest.approx(0.25)


# ---------------------------------------------------------------------------
# osi — the headline correctness checks
# ---------------------------------------------------------------------------

class TestOsi:
    def test_unidirectional_flow_is_zero(self):
        # WSS always points the same way → OSI = 0.
        series = np.tile(np.array([[1.0, 0.0, 0.0]]), (10, 3, 1))  # (10, 3, 3)
        result = osi(series)
        assert result.shape == (3,)
        np.testing.assert_allclose(result, 0.0, atol=1e-12)

    def test_perfectly_reversing_flow_is_half(self):
        # Half the cycle +x, half -x → mean vector zero → OSI = 0.5.
        fwd = np.tile(np.array([[1.0, 0.0, 0.0]]), (5, 4, 1))
        bwd = np.tile(np.array([[-1.0, 0.0, 0.0]]), (5, 4, 1))
        series = np.concatenate([fwd, bwd], axis=0)  # (10, 4, 3)
        np.testing.assert_allclose(osi(series), 0.5, atol=1e-12)

    def test_zero_wss_face_reports_zero_osi(self):
        series = np.zeros((6, 2, 3))
        np.testing.assert_allclose(osi(series), 0.0, atol=1e-12)

    def test_bounded_zero_to_half(self):
        rng = np.random.default_rng(0)
        series = rng.standard_normal((20, 50, 3))
        result = osi(series)
        assert result.min() >= -1e-9
        assert result.max() <= 0.5 + 1e-9

    def test_weights_match_uniform_when_equal(self):
        rng = np.random.default_rng(1)
        series = rng.standard_normal((8, 5, 3))
        w = np.ones(8)
        np.testing.assert_allclose(osi(series, w), osi(series), atol=1e-12)


# ---------------------------------------------------------------------------
# tawss
# ---------------------------------------------------------------------------

class TestTawss:
    def test_constant_magnitude(self):
        # Constant |WSS| = 2 (kinematic) everywhere.
        series = np.tile(np.array([[2.0, 0.0, 0.0]]), (10, 4, 1))
        np.testing.assert_allclose(tawss(series), 2.0)

    def test_is_mean_of_magnitude(self):
        series = np.zeros((2, 1, 3))
        series[0, 0] = [3.0, 0.0, 0.0]   # mag 3
        series[1, 0] = [0.0, 4.0, 0.0]   # mag 4
        np.testing.assert_allclose(tawss(series), 3.5)

    def test_shape(self):
        series = np.ones((5, 7, 3))
        assert tawss(series).shape == (7,)


# ---------------------------------------------------------------------------
# summary_stats — unit conversion and area fractions
# ---------------------------------------------------------------------------

class TestSummaryStats:
    def test_tawss_pa_is_rho_times_kinematic(self):
        tawss_kin = np.array([1.0e-3, 2.0e-3, 3.0e-3])
        osi_field = np.zeros(3)
        areas = np.ones(3)
        stats = summary_stats(tawss_kin, osi_field, areas, rho=RHO)
        # max kinematic 3e-3 → 3e-3 * 1060 Pa
        assert stats["tawss_pa"]["max"] == pytest.approx(3.0e-3 * RHO)
        assert stats["tawss_kinematic"]["max"] == pytest.approx(3.0e-3)

    def test_low_wss_area_fraction(self):
        # Two of four equal-area faces below 0.4 Pa.
        # kinematic = Pa / rho
        tawss_kin = np.array([0.1, 0.2, 1.0, 2.0]) / RHO
        osi_field = np.zeros(4)
        areas = np.ones(4)
        stats = summary_stats(tawss_kin, osi_field, areas, rho=RHO)
        assert stats["area_fraction_tawss_lt_0p4pa"] == pytest.approx(0.5)

    def test_high_osi_area_fraction(self):
        tawss_kin = np.ones(4) * 1e-3
        osi_field = np.array([0.0, 0.1, 0.4, 0.45])  # two faces > 0.3
        areas = np.ones(4)
        stats = summary_stats(tawss_kin, osi_field, areas, rho=RHO)
        assert stats["area_fraction_osi_gt_0p3"] == pytest.approx(0.5)

    def test_area_weighting(self):
        # One big face dominates the weighted mean.
        tawss_kin = np.array([1.0, 0.0])
        osi_field = np.zeros(2)
        areas = np.array([9.0, 1.0])
        stats = summary_stats(tawss_kin, osi_field, areas, rho=1.0)
        assert stats["tawss_pa"]["mean"] == pytest.approx(0.9)

    def test_face_count_and_total_area(self):
        tawss_kin = np.ones(5) * 1e-3
        stats = summary_stats(tawss_kin, np.zeros(5), np.full(5, 2.0), rho=RHO)
        assert stats["n_wall_faces"] == 5
        assert stats["wall_area_m2"] == pytest.approx(10.0)


# ---------------------------------------------------------------------------
# _detect_patches — mode detection from patch_labels.json
# ---------------------------------------------------------------------------

class TestDetectPatches:
    def test_new_mode_detected_when_aneurysm_sac_present(self, tmp_path):
        labels = {"aneurysm_sac": "aneurysm_sac", "parent_vessel": "parent_vessel",
                  "inlet": "inlet", "outlet_0": "outlet"}
        (tmp_path / "patch_labels.json").write_text(json.dumps(labels))
        aneurysm_patch, parent_vessel_patch = _detect_patches(tmp_path)
        assert aneurysm_patch == "aneurysm_sac"
        assert parent_vessel_patch == "parent_vessel"

    def test_legacy_mode_when_wall_present(self, tmp_path):
        labels = {"wall": "wall", "inlet": "inlet", "outlet_0": "outlet"}
        (tmp_path / "patch_labels.json").write_text(json.dumps(labels))
        aneurysm_patch, parent_vessel_patch = _detect_patches(tmp_path)
        assert aneurysm_patch is None
        assert parent_vessel_patch == "wall"

    def test_legacy_mode_when_no_labels_file(self, tmp_path):
        aneurysm_patch, parent_vessel_patch = _detect_patches(tmp_path)
        assert aneurysm_patch is None
        assert parent_vessel_patch == "wall"


# ---------------------------------------------------------------------------
# _read_surface_field_value — CSV parsing
# ---------------------------------------------------------------------------

class TestReadSurfaceFieldValue:
    def _make_dat(self, tmp_path, fo_name, content):
        """Write a surfaceFieldValue.dat file inside a fake postProcessing dir."""
        d = tmp_path / "postProcessing" / fo_name / "0"
        d.mkdir(parents=True)
        (d / "surfaceFieldValue.dat").write_text(textwrap.dedent(content))
        return tmp_path

    def test_scalar_field_parsed(self, tmp_path):
        content = """\
            # Time  p_areaAverage
            0.857   1.234e-3
            0.900   2.345e-3
        """
        self._make_dat(tmp_path, "pressure_mean", content)
        rows = _read_surface_field_value(tmp_path, "pressure_mean", t_start=0.0)
        assert len(rows) == 2
        assert rows[0] == pytest.approx((0.857, 1.234e-3))
        assert rows[1] == pytest.approx((0.900, 2.345e-3))

    def test_vector_field_returns_magnitude(self, tmp_path):
        # Vector (time, vx, vy, vz) → magnitude sqrt(vx²+vy²+vz²)
        content = """\
            # Time  Ux  Uy  Uz
            0.857   3.0  4.0  0.0
        """
        self._make_dat(tmp_path, "neck_flux", content)
        rows = _read_surface_field_value(tmp_path, "neck_flux", t_start=0.0)
        assert len(rows) == 1
        assert rows[0][1] == pytest.approx(5.0)

    def test_parenthesised_vector_is_parsed(self, tmp_path):
        """
        BUG-010 regression. OpenFOAM writes vector values as "(vx vy vz)"; the
        old parser split on whitespace and fed "(3.0" to float(), so every row
        was silently dropped and the metric came out null.
        """
        content = """\
            # Time  (Ux Uy Uz)
            0.857   (3.0 4.0 0.0)
            0.900   (0.0 0.0 5.0)
        """
        self._make_dat(tmp_path, "neck_vel", content)
        rows = _read_surface_field_value(tmp_path, "neck_vel", t_start=0.0)
        assert len(rows) == 2
        assert rows[0][1] == pytest.approx(5.0)
        assert rows[1][1] == pytest.approx(5.0)

    def test_component_reduction_preserves_sign(self, tmp_path):
        """Magnitude would report 2.5 for both rows, hiding the reversal."""
        content = """\
            # Time  (Ux Uy Uz)
            0.857   (-2.5 0.0 0.0)
            0.900   (2.5 0.0 0.0)
        """
        self._make_dat(tmp_path, "flux", content)
        rows = _read_surface_field_value(tmp_path, "flux", t_start=0.0,
                                         reduce="component:0")
        assert rows[0][1] == pytest.approx(-2.5)
        assert rows[1][1] == pytest.approx(2.5)

    def test_malformed_rows_warn_instead_of_silently_dropping(self, tmp_path):
        content = """\
            # Time  value
            0.857   1.0
            0.900   not_a_number
        """
        self._make_dat(tmp_path, "fo", content)
        with pytest.warns(UserWarning, match="unparsable row"):
            rows = _read_surface_field_value(tmp_path, "fo", t_start=0.0)
        assert len(rows) == 1

    def test_unknown_reduction_raises(self, tmp_path):
        content = """\
            0.857   (1.0 2.0 3.0)
        """
        self._make_dat(tmp_path, "fo", content)
        with pytest.raises(ValueError, match="unknown vector reduction"):
            _read_surface_field_value(tmp_path, "fo", t_start=0.0, reduce="bogus")

    def test_comment_lines_skipped(self, tmp_path):
        content = """\
            # header line
            # another comment
            0.857  1.0
        """
        self._make_dat(tmp_path, "fo", content)
        rows = _read_surface_field_value(tmp_path, "fo", t_start=0.0)
        assert len(rows) == 1

    def test_t_start_filter(self, tmp_path):
        content = """\
            0.500  1.0
            0.857  2.0
            1.000  3.0
        """
        self._make_dat(tmp_path, "fo", content)
        rows = _read_surface_field_value(tmp_path, "fo", t_start=0.857)
        assert len(rows) == 2
        assert rows[0][0] == pytest.approx(0.857)

    def test_missing_function_object_returns_empty(self, tmp_path):
        rows = _read_surface_field_value(tmp_path, "nonexistent_fo", t_start=0.0)
        assert rows == []

    def test_merges_multiple_restart_time_dirs(self, tmp_path):
        # A restarted solve writes one time-named subdir per restart. Rows from
        # all of them must be merged in numerical time order (not just the
        # earliest dir), with later dirs winning on overlapping boundary times.
        base = tmp_path / "postProcessing" / "fo"
        (base / "0").mkdir(parents=True)
        (base / "0" / "surfaceFieldValue.dat").write_text(
            textwrap.dedent(
                """\
                0.500  1.0
                0.857  2.0
                """
            )
        )
        # Restart directory; 0.857 overlaps and should be overwritten by 9.0.
        (base / "0.857").mkdir(parents=True)
        (base / "0.857" / "surfaceFieldValue.dat").write_text(
            textwrap.dedent(
                """\
                0.857  9.0
                1.000  3.0
                """
            )
        )
        rows = _read_surface_field_value(tmp_path, "fo", t_start=0.0)
        assert [t for t, _ in rows] == pytest.approx([0.500, 0.857, 1.000])
        assert dict(rows)[0.857] == pytest.approx(9.0)
