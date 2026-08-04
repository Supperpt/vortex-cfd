"""
Tests for vortex_cfd.womersley — pure numpy/scipy, no OpenFOAM required.

The load-bearing property is that each harmonic is area-normalised, so the
area-weighted mean of U(r, t) reproduces U_mean x waveform(t).  Most of these
tests are that identity, checked by radial quadrature.
"""

import numpy as np
import pytest

from vortex_cfd.constants import NU
from vortex_cfd.waveform import T_CYCLE, load_waveform
from vortex_cfd.womersley import (
    _womersley_shape,
    fourier_coefficients,
    womersley_velocities,
    write_boundary_data,
)

RADIUS = 0.003
MEAN_V = 0.4
CENTROID = np.zeros(3)
NORMAL = np.array([0.0, 0.0, 1.0])


@pytest.fixture
def waveform():
    return load_waveform(None)


@pytest.fixture
def radial_faces():
    """
    Faces laid along +x at increasing radius, for radial quadrature:
        <U>(t) = 2 * integral_0^1 U(s, t) * s ds
    Far more accurate than a 2-D grid, and it tests exactly the normalisation.
    """
    n = 400
    s = np.linspace(0.0, 1.0 - 1.0 / n, n)
    faces = np.column_stack([s * RADIUS, np.zeros(n), np.zeros(n)])
    weights = 2.0 * s * (s[1] - s[0])
    return faces, weights


@pytest.fixture
def disc_faces():
    """A centre face plus a ring at r = R/2."""
    theta = np.linspace(0, 2 * np.pi, 8, endpoint=False)
    rim = np.column_stack([
        RADIUS / 2 * np.cos(theta),
        RADIUS / 2 * np.sin(theta),
        np.zeros(8),
    ])
    return np.vstack([np.zeros((1, 3)), rim])


class TestFourierCoefficients:
    def test_shape(self, waveform):
        assert fourier_coefficients(waveform, n_harm=8).shape == (9,)

    def test_dc_term_is_the_waveform_mean(self, waveform):
        # The waveform is normalised to mean 1, so C_0 must be 1.
        c = fourier_coefficients(waveform, n_harm=8)
        assert np.real(c[0]) == pytest.approx(1.0, abs=1e-9)

    def test_is_complex(self, waveform):
        assert fourier_coefficients(waveform, n_harm=4).dtype == complex


class TestWomersleyShape:
    def test_vanishes_at_the_wall(self):
        # G_k(1) = [1 - J0(L)/J0(L)] / N = 0 — the no-slip condition.
        assert abs(_womersley_shape(np.array([1.0]), alpha_k=3.0)[0]) == \
            pytest.approx(0.0, abs=1e-12)

    def test_nonzero_on_the_axis(self):
        assert abs(_womersley_shape(np.array([0.0]), alpha_k=3.0)[0]) > 0.0

    def test_area_mean_is_unity(self):
        """
        The normalisation that makes the harmonics sum to the waveform:
        2 * integral_0^1 G_k(s) s ds == 1.
        """
        n = 2000
        s = np.linspace(0.0, 1.0 - 1.0 / n, n)
        for alpha in (0.5, 3.0, 10.0):
            g = _womersley_shape(s, alpha_k=alpha)
            area_mean = (g * 2.0 * s * (s[1] - s[0])).sum()
            assert area_mean.real == pytest.approx(1.0, abs=0.01)

    def test_low_alpha_tends_to_parabolic(self):
        """As alpha -> 0 the profile relaxes to the quasi-steady Poiseuille shape."""
        s = np.linspace(0, 1, 50)
        g = np.real(_womersley_shape(s, alpha_k=0.05))
        assert g == pytest.approx(2.0 * (1.0 - s ** 2), abs=0.02)


class TestWomersleyVelocities:
    def test_output_shapes(self, waveform, disc_faces):
        t, u = womersley_velocities(disc_faces, CENTROID, NORMAL, RADIUS,
                                    MEAN_V, waveform, NU)
        assert t.shape == (len(waveform) + 1,)
        assert u.shape == (len(waveform) + 1, len(disc_faces), 3)

    def test_velocity_is_parallel_to_the_inward_normal(self, waveform, disc_faces):
        _, u = womersley_velocities(disc_faces, CENTROID, NORMAL, RADIUS,
                                    MEAN_V, waveform, NU)
        assert np.abs(u[:, :, 0]).max() == 0.0
        assert np.abs(u[:, :, 1]).max() == 0.0

    def test_area_mean_reproduces_the_waveform(self, waveform, radial_faces):
        """
        The central identity: area-weighted mean == U_mean * waveform(t), for
        every time step. Deviation is Fourier truncation at 8 harmonics.
        """
        faces, weights = radial_faces
        n = len(waveform)
        _, u = womersley_velocities(faces, CENTROID, NORMAL, RADIUS,
                                    MEAN_V, waveform, NU)
        area_mean = (u[:n, :, 2] * weights[np.newaxis, :]).sum(axis=1)
        expected = waveform[:, 1] * MEAN_V
        assert area_mean == pytest.approx(expected, rel=0.05)

    def test_cycle_mean_matches_mean_velocity(self, waveform, radial_faces):
        faces, weights = radial_faces
        n = len(waveform)
        _, u = womersley_velocities(faces, CENTROID, NORMAL, RADIUS,
                                    MEAN_V, waveform, NU)
        area_mean = (u[:n, :, 2] * weights[np.newaxis, :]).sum(axis=1)
        assert area_mean.mean() == pytest.approx(MEAN_V, rel=1e-3)

    def test_centreline_is_faster_than_the_mean(self, waveform, disc_faces):
        _, u = womersley_velocities(disc_faces, CENTROID, NORMAL, RADIUS,
                                    MEAN_V, waveform, NU)
        assert np.linalg.norm(u[:, 0, :], axis=1).mean() > MEAN_V

    def test_scales_linearly_with_mean_velocity(self, waveform, disc_faces):
        _, u1 = womersley_velocities(disc_faces, CENTROID, NORMAL, RADIUS,
                                     0.3, waveform, NU)
        _, u2 = womersley_velocities(disc_faces, CENTROID, NORMAL, RADIUS,
                                     0.6, waveform, NU)
        np.testing.assert_allclose(u2, 2.0 * u1, rtol=1e-12)

    def test_radius_is_measured_in_the_inlet_plane(self, waveform):
        """
        Faces displaced along the normal are at the same in-plane radius, so
        they must get identical velocities. Catches a missing projection.
        """
        base = np.array([[RADIUS / 2, 0.0, 0.0]])
        shifted = base + np.array([0.0, 0.0, 0.01])   # 10 mm along the normal
        _, u_base = womersley_velocities(base, CENTROID, NORMAL, RADIUS,
                                         MEAN_V, waveform, NU)
        _, u_shift = womersley_velocities(shifted, CENTROID, NORMAL, RADIUS,
                                          MEAN_V, waveform, NU)
        np.testing.assert_allclose(u_base, u_shift, rtol=1e-12)

    def test_rejects_bad_geometry(self, waveform, disc_faces):
        with pytest.raises(ValueError, match="radius must be positive"):
            womersley_velocities(disc_faces, CENTROID, NORMAL, 0.0,
                                 MEAN_V, waveform, NU)
        with pytest.raises(ValueError, match="cycles must be >= 1"):
            womersley_velocities(disc_faces, CENTROID, NORMAL, RADIUS,
                                 MEAN_V, waveform, NU, cycles=0)


class TestMultiCycleCoverage:
    """
    timeVaryingMappedFixedValue does NOT wrap around like flowRateInletVelocity's
    `outOfBounds repeat`: past the last supplied time it holds the final value.
    The prototype branch wrote one cycle only, so cycles > 1 would have frozen
    the inlet after the first cycle. These tests pin the fix.
    """

    @pytest.mark.parametrize("cycles", [1, 2, 3, 5])
    def test_data_spans_the_whole_run(self, waveform, disc_faces, cycles):
        t, _ = womersley_velocities(disc_faces, CENTROID, NORMAL, RADIUS,
                                    MEAN_V, waveform, NU, cycles=cycles)
        assert t[-1] >= cycles * T_CYCLE - 1e-12

    @pytest.mark.parametrize("cycles", [1, 2, 3])
    def test_times_are_strictly_increasing(self, waveform, disc_faces, cycles):
        t, _ = womersley_velocities(disc_faces, CENTROID, NORMAL, RADIUS,
                                    MEAN_V, waveform, NU, cycles=cycles)
        assert np.all(np.diff(t) > 0)

    def test_sample_count_grows_with_cycles(self, waveform, disc_faces):
        n = len(waveform)
        for cycles in (1, 2, 3):
            t, u = womersley_velocities(disc_faces, CENTROID, NORMAL, RADIUS,
                                        MEAN_V, waveform, NU, cycles=cycles)
            assert len(t) == cycles * n + 1
            assert len(u) == cycles * n + 1

    def test_profile_repeats_each_cycle(self, waveform, disc_faces):
        n = len(waveform)
        _, u = womersley_velocities(disc_faces, CENTROID, NORMAL, RADIUS,
                                    MEAN_V, waveform, NU, cycles=3)
        np.testing.assert_allclose(u[:n], u[n:2 * n], rtol=1e-12)
        np.testing.assert_allclose(u[:n], u[2 * n:3 * n], rtol=1e-12)

    def test_terminal_sample_closes_the_cycle(self, waveform, disc_faces):
        """The appended end-point must equal t=0 so the seam is continuous."""
        _, u = womersley_velocities(disc_faces, CENTROID, NORMAL, RADIUS,
                                    MEAN_V, waveform, NU, cycles=2)
        np.testing.assert_allclose(u[-1], u[0], rtol=1e-12)


class TestWriteBoundaryData:
    def _write(self, tmp_path, waveform, faces, cycles=1):
        t, u = womersley_velocities(faces, CENTROID, NORMAL, RADIUS,
                                    MEAN_V, waveform, NU, cycles=cycles)
        write_boundary_data(tmp_path, "inlet", t, faces, u)
        return t, u, tmp_path / "constant" / "boundaryData" / "inlet"

    def test_writes_points_and_one_u_per_time(self, tmp_path, waveform, disc_faces):
        t, _, bd = self._write(tmp_path, waveform, disc_faces)
        assert (bd / "points").exists()
        assert len(list(bd.glob("*/U"))) == len(t)

    def test_points_count_matches_faces(self, tmp_path, waveform, disc_faces):
        _, _, bd = self._write(tmp_path, waveform, disc_faces)
        lines = (bd / "points").read_text().splitlines()
        assert str(len(disc_faces)) in lines

    def test_u_files_have_one_vector_per_face(self, tmp_path, waveform, disc_faces):
        _, _, bd = self._write(tmp_path, waveform, disc_faces)
        text = next(iter(bd.glob("*/U"))).read_text()
        assert text.count("(") == len(disc_faces) + 1   # +1 for the list opener

    def test_multicycle_directories_reach_end_time(self, tmp_path, waveform, disc_faces):
        _, _, bd = self._write(tmp_path, waveform, disc_faces, cycles=3)
        written = sorted(float(p.parent.name) for p in bd.glob("*/U"))
        assert max(written) >= 3 * T_CYCLE - 1e-9

    def test_directory_names_parse_as_times(self, tmp_path, waveform, disc_faces):
        t, _, bd = self._write(tmp_path, waveform, disc_faces, cycles=2)
        written = sorted(float(p.parent.name) for p in bd.glob("*/U"))
        # Names are formatted to 5 dp, so compare at that resolution.
        np.testing.assert_allclose(written, np.round(t, 5), atol=1e-5)

    def test_mismatched_lengths_raise(self, tmp_path, disc_faces):
        with pytest.raises(ValueError, match="disagree"):
            write_boundary_data(tmp_path, "inlet", np.array([0.0, 1.0]),
                                disc_faces, np.zeros((1, len(disc_faces), 3)))
