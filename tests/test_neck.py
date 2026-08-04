"""
Tests for vortex_cfd.neck — pure numpy/pyvista, no OpenFOAM required.

Two separable concerns, tested separately:
  * geometry   — fitting the orifice plane to the sac's open boundary loop
  * quadrature — integrating flux over the clipped disc, against closed forms
"""

import json

import numpy as np
import pytest
import pyvista as pv

from vortex_cfd.neck import (
    NeckGeometryError,
    NeckPlane,
    disc_sample,
    extract_neck_orifice,
    flux_metrics,
    load_neck_plane,
    resolve_neck_plane,
)

from tests.conftest import _open_sac_stl, _sphere_stl

SAC_R = 0.004
CLIP_Z = 0.002
# A sphere of radius R cut at height z has a rim of radius sqrt(R^2 - z^2).
EXACT_LOOP_R = np.sqrt(SAC_R ** 2 - CLIP_Z ** 2)


def _plane(radius, origin=(0, 0, 0), normal=(0, 0, 1)):
    """A NeckPlane built directly, bypassing STL extraction."""
    return NeckPlane(
        origin=np.asarray(origin, dtype=float),
        normal=np.asarray(normal, dtype=float),
        radius=radius,
        r_min=radius,
        r_max=radius,
        r_eff=radius,
        loop_area_m2=np.pi * radius ** 2,
        planarity=0.0,
        n_loop_points=48,
        area_fraction_positive_side=1.0,
    )


def _box_grid(n, length=0.01, u_field=None):
    """Uniform cubic grid centred on the origin, carrying cell-data U."""
    edges = np.linspace(-length / 2, length / 2, n + 1)
    grid = pv.RectilinearGrid(edges, edges, edges).cast_to_unstructured_grid()
    if u_field is None:
        u_field = np.tile([0.0, 0.0, 1.0], (grid.n_cells, 1))
    grid.cell_data["U"] = u_field
    return grid


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

class TestExtractNeckOrifice:
    def test_origin_and_radius(self, tmp_path):
        p = _open_sac_stl(tmp_path / "sac.stl", radius=SAC_R, clip_z=CLIP_Z)
        plane = extract_neck_orifice(p)
        assert plane.origin == pytest.approx([0.0, 0.0, CLIP_Z], abs=1e-5)
        # r_eff is deliberately a slight under-estimate (polygon inscribed in
        # the true circle); 2 % covers it with margin.
        assert plane.r_eff == pytest.approx(EXACT_LOOP_R, rel=0.02)
        assert plane.radius == plane.r_eff
        assert plane.r_min <= plane.r_eff <= plane.r_max

    def test_loop_is_planar(self, tmp_path):
        p = _open_sac_stl(tmp_path / "sac.stl", radius=SAC_R, clip_z=CLIP_Z)
        plane = extract_neck_orifice(p)
        assert plane.planarity < 1e-4
        assert plane.n_loop_points > 10

    def test_normal_is_unit_and_points_into_sac(self, tmp_path):
        p = _open_sac_stl(tmp_path / "sac.stl", radius=SAC_R, clip_z=CLIP_Z)
        plane = extract_neck_orifice(p)
        assert np.linalg.norm(plane.normal) == pytest.approx(1.0)
        assert plane.normal == pytest.approx([0.0, 0.0, 1.0], abs=1e-6)

    def test_normal_flips_when_sac_is_on_the_other_side(self, tmp_path):
        """
        The single most important test in this module.  The normal's sign
        decides whether max(U.n, 0) measures inflow or outflow; a flipped
        normal yields a plausible-looking number for the wrong quantity.
        It must be derived from the sac geometry, never from neck_plane.json.
        """
        p = _open_sac_stl(tmp_path / "sac.stl", radius=SAC_R, clip_z=CLIP_Z,
                          keep_above=False)
        plane = extract_neck_orifice(p)
        assert plane.normal == pytest.approx([0.0, 0.0, -1.0], abs=1e-6)
        assert plane.area_fraction_positive_side > 0.95

    def test_oblique_clip_is_fitted_correctly(self, tmp_path):
        n = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0)
        p = _open_sac_stl(tmp_path / "sac.stl", radius=SAC_R, clip_z=CLIP_Z,
                          normal=tuple(n))
        plane = extract_neck_orifice(p)
        assert plane.normal == pytest.approx(n, abs=1e-3)
        assert plane.r_eff == pytest.approx(EXACT_LOOP_R, rel=0.02)

    def test_closed_sac_raises(self, tmp_path):
        p = _sphere_stl(tmp_path / "closed.stl", radius=SAC_R)
        with pytest.raises(NeckGeometryError, match="closed surface"):
            extract_neck_orifice(p)

    def test_two_loops_raises(self, tmp_path):
        mesh = pv.Sphere(radius=SAC_R, theta_resolution=24, phi_resolution=24)
        mesh = mesh.clip(normal=(0, 0, 1), origin=(0, 0, CLIP_Z), invert=True)
        mesh = mesh.clip(normal=(0, 0, 1), origin=(0, 0, -CLIP_Z), invert=False)
        path = tmp_path / "band.stl"
        mesh.extract_surface(algorithm="dataset_surface").save(str(path), binary=False)
        with pytest.raises(NeckGeometryError, match="2 open boundary loops"):
            extract_neck_orifice(path)


class TestResolveNeckPlane:
    def _write_json(self, path, origin, normal):
        path.write_text(json.dumps({"origin": list(origin), "normal": list(normal)}))
        return path

    def test_cross_check_agrees(self, tmp_path):
        sac = _open_sac_stl(tmp_path / "sac.stl", radius=SAC_R, clip_z=CLIP_Z)
        js = self._write_json(tmp_path / "neck_plane.json", [0, 0, CLIP_Z], [0, 0, 1])
        out = resolve_neck_plane(sac, js)
        assert out["cross_check"]["agrees"] is True
        assert out["cross_check"]["json_normal_angle_deg"] == pytest.approx(0.0, abs=1e-6)
        assert out["cross_check"]["json_normal_sign_agrees"] is True

    def test_reversed_json_normal_does_not_flip_the_fit(self, tmp_path):
        """neck_plane.json is a cross-check only — it must not steer the fit."""
        sac = _open_sac_stl(tmp_path / "sac.stl", radius=SAC_R, clip_z=CLIP_Z)
        js = self._write_json(tmp_path / "neck_plane.json", [0, 0, CLIP_Z], [0, 0, -1])
        out = resolve_neck_plane(sac, js)
        assert out["normal"] == pytest.approx([0.0, 0.0, 1.0], abs=1e-6)
        assert out["cross_check"]["json_normal_sign_agrees"] is False

    def test_non_unit_json_normal_is_reported(self, tmp_path):
        sac = _open_sac_stl(tmp_path / "sac.stl", radius=SAC_R, clip_z=CLIP_Z)
        js = self._write_json(tmp_path / "neck_plane.json", [0, 0, CLIP_Z], [0, 0, 7.5])
        out = resolve_neck_plane(sac, js)
        assert out["cross_check"]["json_normal_was_unit"] is False
        assert out["cross_check"]["json_normal_angle_deg"] == pytest.approx(0.0, abs=1e-6)

    def test_missing_json_still_resolves(self, tmp_path):
        sac = _open_sac_stl(tmp_path / "sac.stl", radius=SAC_R, clip_z=CLIP_Z)
        out = resolve_neck_plane(sac, None)
        assert out["cross_check"] is None
        assert out["radius_m"] > 0

    def test_roundtrip_through_json(self, tmp_path):
        sac = _open_sac_stl(tmp_path / "sac.stl", radius=SAC_R, clip_z=CLIP_Z)
        out = resolve_neck_plane(sac, None)
        (tmp_path / "neck_plane_resolved.json").write_text(json.dumps(out))
        plane = load_neck_plane(tmp_path)
        assert plane is not None
        assert plane.radius == pytest.approx(out["radius_m"])
        assert plane.normal == pytest.approx(out["normal"])

    def test_load_returns_none_when_absent(self, tmp_path):
        assert load_neck_plane(tmp_path) is None


# ---------------------------------------------------------------------------
# Pure metrics — closed forms, no mesh
# ---------------------------------------------------------------------------

class TestFluxMetrics:
    N = np.array([0.0, 0.0, 1.0])

    def test_balanced_flow_gives_zero_net_but_positive_inflow(self):
        areas = np.array([1.0, 1.0])
        vecs = np.array([[0, 0, 2.0], [0, 0, -2.0]])
        m = flux_metrics(areas, vecs, self.N)
        assert m["inflow_rate_m3s"] == pytest.approx(2.0)
        assert m["outflow_rate_m3s"] == pytest.approx(-2.0)
        assert m["net_flux_m3s"] == pytest.approx(0.0, abs=1e-15)

    def test_all_outflow_gives_zero_inflow(self):
        """
        BUG-011 regression: the old code took max() of a signed flux series, so
        an all-negative series returned the least-negative value as the 'peak'.
        A positive-part integral cannot do that.
        """
        areas = np.array([1.0, 1.0, 1.0])
        vecs = np.array([[0, 0, -1.0], [0, 0, -3.0], [0, 0, -2.0]])
        m = flux_metrics(areas, vecs, self.N)
        assert m["inflow_rate_m3s"] == 0.0
        assert m["net_flux_m3s"] < 0
        assert m["peak_inflow_velocity_ms"] == 0.0

    def test_peak_velocity_counts_outward_cells_but_peak_inflow_does_not(self):
        areas = np.array([1.0, 1.0])
        vecs = np.array([[0, 0, 1.0], [0, 0, -9.0]])
        m = flux_metrics(areas, vecs, self.N)
        assert m["peak_velocity_ms"] == pytest.approx(9.0)
        assert m["peak_inflow_velocity_ms"] == pytest.approx(1.0)

    def test_tangential_flow_carries_no_flux(self):
        areas = np.array([1.0])
        vecs = np.array([[5.0, 5.0, 0.0]])
        m = flux_metrics(areas, vecs, self.N)
        assert m["net_flux_m3s"] == pytest.approx(0.0, abs=1e-15)
        assert m["inflow_rate_m3s"] == pytest.approx(0.0, abs=1e-15)
        assert m["peak_velocity_ms"] > 0

    def test_normal_is_normalised_internally(self):
        areas = np.array([1.0])
        vecs = np.array([[0, 0, 2.0]])
        m = flux_metrics(areas, vecs, np.array([0.0, 0.0, 50.0]))
        assert m["net_flux_m3s"] == pytest.approx(2.0)

    def test_empty_disc_is_zero_and_warns(self):
        with pytest.warns(UserWarning, match="empty neck disc"):
            m = flux_metrics(np.zeros(0), np.zeros((0, 3)), self.N)
        assert m["inflow_rate_m3s"] == 0.0
        assert m["n_cells"] == 0


# ---------------------------------------------------------------------------
# Quadrature over the clipped disc — against analytic flows
# ---------------------------------------------------------------------------

class TestDiscSample:
    R = 0.002

    def test_uniform_flow_is_exact_against_measured_area(self):
        """Discretisation cancels: Q must equal (measured area) x speed."""
        grid = _box_grid(40)
        areas, vecs = disc_sample(grid, _plane(self.R))
        m = flux_metrics(areas, vecs, np.array([0.0, 0.0, 1.0]))
        assert m["net_flux_m3s"] == pytest.approx(m["disc_area_m2"], rel=1e-12)
        assert m["inflow_rate_m3s"] == pytest.approx(m["disc_area_m2"], rel=1e-12)

    def test_disc_area_converges_to_pi_r_squared(self):
        grid = _box_grid(40)
        areas, _ = disc_sample(grid, _plane(self.R))
        assert areas.sum() == pytest.approx(np.pi * self.R ** 2, rel=0.02)

    def test_parabolic_profile_recovers_mean_flow_rate(self):
        grid = _box_grid(40)
        c = grid.cell_centers().points
        r = np.hypot(c[:, 0], c[:, 1])
        u_mean = 0.3
        uz = 2 * u_mean * np.clip(1 - (r / self.R) ** 2, 0, None)
        grid.cell_data["U"] = np.column_stack([np.zeros_like(uz), np.zeros_like(uz), uz])

        areas, vecs = disc_sample(grid, _plane(self.R))
        m = flux_metrics(areas, vecs, np.array([0.0, 0.0, 1.0]))
        assert m["inflow_rate_m3s"] == pytest.approx(u_mean * np.pi * self.R ** 2, rel=0.02)
        # Cell-centre sampling never lands exactly on the axis, so the peak is a
        # slight under-estimate — systematic, and it shrinks with refinement.
        assert m["peak_velocity_ms"] == pytest.approx(2 * u_mean, rel=0.05)

    def test_reversing_flow_nets_to_zero_exactly(self):
        grid = _box_grid(40)
        c = grid.cell_centers().points
        uz = np.where(c[:, 0] > 0, 1.0, -1.0)
        grid.cell_data["U"] = np.column_stack([np.zeros_like(uz), np.zeros_like(uz), uz])

        areas, vecs = disc_sample(grid, _plane(self.R))
        m = flux_metrics(areas, vecs, np.array([0.0, 0.0, 1.0]))
        assert m["net_flux_m3s"] == pytest.approx(0.0, abs=1e-18)
        assert m["inflow_rate_m3s"] == pytest.approx(m["disc_area_m2"] / 2, rel=1e-12)

    def test_oblique_plane_projects_radius_in_plane(self):
        """
        A missing in-plane projection of the radius silently mis-sizes the disc
        on any non-axis-aligned plane.  The flux-per-area invariant is exact
        regardless of discretisation, so it isolates that bug.
        """
        n = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0)
        grid = _box_grid(60, u_field=None)
        grid.cell_data["U"] = np.tile(n, (grid.n_cells, 1))
        areas, vecs = disc_sample(grid, _plane(self.R, normal=tuple(n)))
        m = flux_metrics(areas, vecs, n)
        assert m["n_cells"] > 50
        assert m["net_flux_m3s"] == pytest.approx(m["disc_area_m2"], rel=1e-12)

    def test_plane_outside_mesh_warns_and_returns_empty(self):
        grid = _box_grid(20)
        with pytest.warns(UserWarning, match="does not intersect"):
            areas, vecs = disc_sample(grid, _plane(self.R, origin=(0, 0, 1.0)))
        assert areas.size == 0
        assert vecs.shape == (0, 3)

    def test_subcell_radius_warns(self):
        grid = _box_grid(20)
        with pytest.warns(UserWarning):
            areas, _ = disc_sample(grid, _plane(1e-7))
        assert areas.size < 20

    def test_missing_field_raises_clearly(self):
        grid = _box_grid(20)
        with pytest.raises(KeyError, match="not found"):
            disc_sample(grid, _plane(self.R), field="NoSuchField")
