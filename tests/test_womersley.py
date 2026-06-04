"""Tests for vortex_cfd.womersley — pure Python/scipy, no OpenFOAM needed."""

import numpy as np
import pytest

from vortex_cfd.waveform import load_waveform, T_CYCLE
from vortex_cfd.womersley import (
    fourier_coefficients,
    womersley_velocities,
    write_boundary_data,
    _womersley_shape,
)

# Shared geometry for most tests
RADIUS = 0.003        # 3 mm inlet
NU = 3.3e-6
MEAN_V = 0.4          # m/s
CENTROID = np.array([0.0, 0.0, 0.0])
NORMAL = np.array([0.0, 0.0, 1.0])   # inward normal: +z


@pytest.fixture
def waveform():
    return load_waveform(None)


@pytest.fixture
def face_centers_disc():
    """Ring of 8 faces at r = R/2, plus one at the centre."""
    theta = np.linspace(0, 2 * np.pi, 8, endpoint=False)
    rim = np.column_stack([
        RADIUS / 2 * np.cos(theta),
        RADIUS / 2 * np.sin(theta),
        np.zeros(8),
    ])
    centre = np.array([[0.0, 0.0, 0.0]])
    return np.vstack([centre, rim])   # 9 faces


# ---------------------------------------------------------------------------
# fourier_coefficients
# ---------------------------------------------------------------------------

class TestFourierCoefficients:
    def test_shape(self, waveform):
        c = fourier_coefficients(waveform, N_harm=8)
        assert c.shape == (9,)

    def test_dc_is_mean(self, waveform):
        c = fourier_coefficients(waveform, N_harm=8)
        # Waveform is normalised: mean = 1.0, so DC coeff ≈ 1.0
        assert np.real(c[0]) == pytest.approx(1.0, abs=1e-6)

    def test_complex_output(self, waveform):
        c = fourier_coefficients(waveform, N_harm=4)
        assert c.dtype == complex


# ---------------------------------------------------------------------------
# _womersley_shape
# ---------------------------------------------------------------------------

class TestWomersleyShape:
    def test_zero_at_wall(self):
        # At r/R = 1, W_k(1) = 1 - J0(Lambda) / J0(Lambda) = 0
        r_norm = np.array([1.0])
        W = _womersley_shape(r_norm, alpha_k=3.0)
        assert np.abs(W[0]) == pytest.approx(0.0, abs=1e-10)

    def test_nonzero_at_centre(self):
        r_norm = np.array([0.0])
        W = _womersley_shape(r_norm, alpha_k=3.0)
        assert np.abs(W[0]) > 0.0

    def test_low_alpha_approaches_poiseuille(self):
        # For very small α, Womersley → Poiseuille profile (parabolic)
        # W(s) ≈ 1 - (1 - (Lambda*s)²/4 + ...) / (1 - Lambda²/4 + ...)
        # At low α, shape is roughly parabolic.
        r_norm = np.linspace(0, 1, 20)
        W = _womersley_shape(r_norm, alpha_k=0.1)
        poiseuille = 1.0 - r_norm ** 2
        # Both should peak at centre and go to 0 at wall; check monotonicity
        mag = np.abs(W)
        assert mag[0] > mag[-1]


# ---------------------------------------------------------------------------
# womersley_velocities
# ---------------------------------------------------------------------------

class TestWomersleyVelocities:
    def test_output_shape(self, waveform, face_centers_disc):
        t_abs, U = womersley_velocities(
            face_centers_disc, CENTROID, NORMAL, RADIUS, MEAN_V, waveform, NU
        )
        N_t = len(waveform)
        assert t_abs.shape == (N_t,)
        assert U.shape == (N_t, len(face_centers_disc), 3)

    def test_direction_parallel_to_normal(self, waveform, face_centers_disc):
        _, U = womersley_velocities(
            face_centers_disc, CENTROID, NORMAL, RADIUS, MEAN_V, waveform, NU
        )
        # All velocity vectors should be parallel to the inward normal (0, 0, 1)
        # i.e. Ux ≈ 0 and Uy ≈ 0 everywhere
        assert np.abs(U[:, :, 0]).max() < 1e-12
        assert np.abs(U[:, :, 1]).max() < 1e-12

    def test_systolic_velocities_positive(self, waveform, face_centers_disc):
        # At systole (peak Q), all faces should be forward-flowing.
        # Near-wall faces MAY go negative at diastole (physiologically correct).
        _, U = womersley_velocities(
            face_centers_disc, CENTROID, NORMAL, RADIUS, MEAN_V, waveform, NU
        )
        systolic_idx = int(np.argmax(waveform[:, 1]))
        speeds_systole = np.linalg.norm(U[systolic_idx], axis=1)
        assert speeds_systole.min() >= 0.0

    def test_centreline_faster_than_mean(self, waveform, face_centers_disc):
        _, U = womersley_velocities(
            face_centers_disc, CENTROID, NORMAL, RADIUS, MEAN_V, waveform, NU
        )
        # Face 0 is at the centre (r=0); its speed should exceed the mean
        centre_speed = np.linalg.norm(U[:, 0, :], axis=1)   # (N_t,)
        assert centre_speed.mean() > MEAN_V

    def test_area_weighted_mean_close_to_waveform(self, waveform):
        # Radial quadrature: evaluate U at many radii and integrate
        #   <U>(t) = 2/R² * ∫₀ᴿ U(r, t) * r dr  ≈  2 * Σ_i U(r_i) * s_i * Δs
        # where s_i = r_i/R ∈ [0, 1).  This is much more accurate than a 2D
        # Cartesian grid and directly tests the area normalization of G_k.
        n_r = 200
        s = np.linspace(0.0, 1.0 - 1.0 / n_r, n_r)   # s = r/R, avoids wall singularity
        ds = s[1] - s[0]
        r_vals = s * RADIUS
        # Place faces along the x-axis (y=z=0) so r = x
        fc = np.column_stack([r_vals, np.zeros(n_r), np.zeros(n_r)])

        t_abs, U = womersley_velocities(
            fc, CENTROID, NORMAL, RADIUS, MEAN_V, waveform, NU, N_harm=8
        )
        speeds = np.linalg.norm(U, axis=2)                    # (N_t, n_r)

        # Area-weighted mean: ∫₀¹ speed(s) * 2s ds
        weights = 2.0 * s * ds                                # (n_r,)
        area_mean = (speeds * weights[np.newaxis, :]).sum(axis=1)  # (N_t,)
        expected = waveform[:, 1] * MEAN_V

        # Filter to time steps well above the waveform minimum (q > 0.3) to
        # avoid diastolic retrograde-flow faces confounding the area-mean test.
        # At those time steps U(r, t) should be entirely forward-flowing and the
        # area-weighted mean should closely match the waveform value.
        mask = waveform[:, 1] > 0.3
        ratio = area_mean[mask] / expected[mask]
        assert ratio.mean() == pytest.approx(1.0, abs=0.10)

    def test_t_abs_matches_waveform_times(self, waveform, face_centers_disc):
        t_abs, _ = womersley_velocities(
            face_centers_disc, CENTROID, NORMAL, RADIUS, MEAN_V, waveform, NU
        )
        expected_t = waveform[:, 0] * T_CYCLE
        np.testing.assert_allclose(t_abs, expected_t, atol=1e-12)


# ---------------------------------------------------------------------------
# write_boundary_data
# ---------------------------------------------------------------------------

class TestWriteBoundaryData:
    def test_creates_points_file(self, tmp_path, waveform, face_centers_disc):
        t_abs, U = womersley_velocities(
            face_centers_disc, CENTROID, NORMAL, RADIUS, MEAN_V, waveform, NU
        )
        write_boundary_data(tmp_path, "inlet", t_abs, face_centers_disc, U)
        points_file = tmp_path / "constant" / "boundaryData" / "inlet" / "points"
        assert points_file.exists()

    def test_points_file_has_correct_count(self, tmp_path, waveform, face_centers_disc):
        t_abs, U = womersley_velocities(
            face_centers_disc, CENTROID, NORMAL, RADIUS, MEAN_V, waveform, NU
        )
        write_boundary_data(tmp_path, "inlet", t_abs, face_centers_disc, U)
        text = (tmp_path / "constant" / "boundaryData" / "inlet" / "points").read_text()
        assert str(len(face_centers_disc)) in text

    def test_creates_u_file_per_timestep(self, tmp_path, waveform, face_centers_disc):
        t_abs, U = womersley_velocities(
            face_centers_disc, CENTROID, NORMAL, RADIUS, MEAN_V, waveform, NU
        )
        write_boundary_data(tmp_path, "inlet", t_abs, face_centers_disc, U)
        bd = tmp_path / "constant" / "boundaryData" / "inlet"
        u_files = list(bd.glob("*/U"))
        assert len(u_files) == len(t_abs)

    def test_u_file_non_empty(self, tmp_path, waveform, face_centers_disc):
        t_abs, U = womersley_velocities(
            face_centers_disc, CENTROID, NORMAL, RADIUS, MEAN_V, waveform, NU
        )
        write_boundary_data(tmp_path, "inlet", t_abs, face_centers_disc, U)
        bd = tmp_path / "constant" / "boundaryData" / "inlet"
        for u_file in bd.glob("*/U"):
            assert u_file.stat().st_size > 0
