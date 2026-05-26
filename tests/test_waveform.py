"""Tests for vortex_cfd.waveform."""

import csv
import textwrap

import numpy as np
import pytest

from vortex_cfd.waveform import N_POINTS, T_CYCLE, load_waveform


class TestDefaultWaveform:
    def test_returns_ndarray(self):
        wf = load_waveform(None)
        assert isinstance(wf, np.ndarray)

    def test_shape(self):
        wf = load_waveform(None)
        assert wf.shape == (N_POINTS, 2)

    def test_time_starts_at_zero(self):
        wf = load_waveform(None)
        assert wf[0, 0] == pytest.approx(0.0)

    def test_time_does_not_reach_one(self):
        # endpoint=False — the table covers [0, 1) so OpenFOAM's repeat is seamless
        wf = load_waveform(None)
        assert wf[-1, 0] < 1.0

    def test_flow_mean_is_one(self):
        wf = load_waveform(None)
        assert np.mean(wf[:, 1]) == pytest.approx(1.0, rel=1e-6)

    def test_flow_all_positive(self):
        # ICA flow should be antegrade throughout the cycle (clipped to > 0)
        wf = load_waveform(None)
        assert np.all(wf[:, 1] > 0)

    def test_systolic_peak_near_30_percent(self):
        wf = load_waveform(None)
        peak_t = wf[np.argmax(wf[:, 1]), 0]
        assert 0.20 <= peak_t <= 0.45, (
            f"Systolic peak at t_norm={peak_t:.2f} — expected in [0.20, 0.45]"
        )

    def test_systolic_peak_exceeds_mean(self):
        wf = load_waveform(None)
        assert np.max(wf[:, 1]) > 1.5, "Systolic peak should be well above the mean"


class TestCSVWaveform:
    def test_loads_csv(self, tmp_path):
        p = tmp_path / "wf.csv"
        p.write_text("0.0,1.0\n0.5,3.0\n1.0,1.0\n")
        wf = load_waveform(str(p))
        assert wf.shape == (3, 2)

    def test_normalises_to_mean_one(self, tmp_path):
        p = tmp_path / "wf.csv"
        p.write_text("0.0,2.0\n0.5,4.0\n1.0,6.0\n")
        wf = load_waveform(str(p))
        assert np.mean(wf[:, 1]) == pytest.approx(1.0, rel=1e-6)

    def test_skips_comment_lines(self, tmp_path):
        p = tmp_path / "wf.csv"
        p.write_text("# time_norm, flow_norm\n0.0,1.0\n0.5,2.0\n")
        wf = load_waveform(str(p))
        assert wf.shape[0] == 2

    def test_skips_blank_lines(self, tmp_path):
        p = tmp_path / "wf.csv"
        p.write_text("0.0,1.0\n\n0.5,2.0\n")
        wf = load_waveform(str(p))
        assert wf.shape[0] == 2

    def test_raises_if_too_short(self, tmp_path):
        p = tmp_path / "wf.csv"
        p.write_text("0.0,1.0\n")
        with pytest.raises(ValueError, match="at least 2"):
            load_waveform(str(p))

    def test_time_column_preserved(self, tmp_path):
        p = tmp_path / "wf.csv"
        p.write_text("0.0,1.0\n0.25,2.0\n0.75,1.5\n")
        wf = load_waveform(str(p))
        np.testing.assert_array_almost_equal(wf[:, 0], [0.0, 0.25, 0.75])

    def test_constant_flow_normalises_to_one(self, tmp_path):
        p = tmp_path / "wf.csv"
        p.write_text("0.0,5.0\n0.5,5.0\n1.0,5.0\n")
        wf = load_waveform(str(p))
        np.testing.assert_array_almost_equal(wf[:, 1], [1.0, 1.0, 1.0])
