"""
Neck-orifice geometry and inflow metrics for the aneurysm sac.

**Why this is not an OpenFOAM function object.**  The clinical quantity is the
*neck inflow rate*: the volume of blood entering the sac per second, i.e. the
integral of only the **inward** part of U.n over the neck orifice.  A
``surfaceFieldValue`` FO can compute ``areaNormalIntegrate`` (the *net* flux) but
cannot express a positive-part integral, and for a sealed sac the net cyclic flux
is ~0 by conservation.  So the metrics are computed here in Python by slicing the
saved U volume snapshots.  This also fixed CAVEAT-012: the old FO used an
*infinite* sampling plane that integrated the whole parent-vessel cross-section.

**Where the plane comes from.**  Entirely from ``aneurysm_sac.stl``.  VORTEX clips
the sac at the neck, so the sac's open boundary loop *is* the orifice: fitting a
plane to that loop gives origin, normal, and extent that are mutually consistent
by construction.  ``neck_plane.json`` is used only as a cross-check, never for
correctness — its normal is not guaranteed to be unit-length and its sign
convention is undocumented, and a flipped normal would silently turn the reported
inflow into the *outflow* (a plausible-looking number that is simply the wrong
metric).  The normal's sign is therefore derived geometrically, by orienting it
towards the sac.

Layout mirrors ``postprocess``: pure functions (``flux_metrics``) carry the unit
tests, the pyvista layer (``disc_sample``, ``read_neck_series``) does the I/O.
Time-averaging deliberately lives in ``postprocess`` so this stays a leaf module.
"""
from __future__ import annotations

import json
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pyvista as pv

from .scaling import read_stl

# A loop this far from planar is suspicious but not fatal — VMTK clips can be
# mildly non-planar.  Ratio of smallest to largest SVD singular value.
PLANARITY_WARN = 0.02
# Below this many cells the disc quadrature is meaningless.
MIN_DISC_CELLS = 20
# Fraction of sac area that must lie on the +normal side for the orientation to
# be unambiguous (a clipped dome is essentially all on one side).
MIN_AREA_FRACTION = 0.95

RESOLVED_FILENAME = "neck_plane_resolved.json"


class NeckGeometryError(ValueError):
    """The neck orifice could not be determined from the aneurysm sac STL."""


@dataclass(frozen=True, eq=False)
class NeckPlane:
    """
    The neck orifice, fitted to the open boundary loop of the sac STL.

    origin  : loop centroid [m]
    normal  : unit normal, oriented *into* the sac
    radius  : the disc radius used for sampling — ``r_eff`` by default
    r_eff   : area-equivalent radius sqrt(A_loop / pi); preferred over r_max
              because over-inclusion re-creates CAVEAT-012 (the disc reaches
              into the parent vessel) while under-inclusion is a mild area bias
    """
    origin: np.ndarray
    normal: np.ndarray
    radius: float
    r_min: float
    r_max: float
    r_eff: float
    loop_area_m2: float
    planarity: float
    n_loop_points: int
    area_fraction_positive_side: float
    source: str = "aneurysm_sac.stl open boundary loop"

    def to_dict(self) -> dict:
        return {
            "origin": [float(v) for v in self.origin],
            "normal": [float(v) for v in self.normal],
            "radius_m": float(self.radius),
            "r_min": float(self.r_min),
            "r_max": float(self.r_max),
            "r_eff": float(self.r_eff),
            "loop_area_m2": float(self.loop_area_m2),
            "planarity": float(self.planarity),
            "n_loop_points": int(self.n_loop_points),
            "area_fraction_positive_side": float(self.area_fraction_positive_side),
            "source": self.source,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "NeckPlane":
        return cls(
            origin=np.asarray(d["origin"], dtype=float),
            normal=np.asarray(d["normal"], dtype=float),
            radius=float(d["radius_m"]),
            r_min=float(d.get("r_min", d["radius_m"])),
            r_max=float(d.get("r_max", d["radius_m"])),
            r_eff=float(d.get("r_eff", d["radius_m"])),
            loop_area_m2=float(d.get("loop_area_m2", 0.0)),
            planarity=float(d.get("planarity", 0.0)),
            n_loop_points=int(d.get("n_loop_points", 0)),
            area_fraction_positive_side=float(d.get("area_fraction_positive_side", 1.0)),
            source=d.get("source", "unknown"),
        )


# ---------------------------------------------------------------------------
# Geometry — fitting the orifice to the sac's open boundary loop
# ---------------------------------------------------------------------------

def _area_weighted_centroid(mesh) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (centroid, cell_centres, cell_areas) for a surface mesh."""
    sized = mesh.compute_cell_sizes(length=False, area=True, volume=False)
    areas = np.asarray(sized.cell_data["Area"], dtype=float)
    centres = np.asarray(mesh.cell_centers().points, dtype=float)
    total = areas.sum()
    if total <= 0:
        raise NeckGeometryError(
            "aneurysm sac STL has zero total surface area — the file is empty "
            "or degenerate."
        )
    return (areas[:, None] * centres).sum(axis=0) / total, centres, areas


def _loop_area_and_radii(
    points: np.ndarray, origin: np.ndarray, e1: np.ndarray, e2: np.ndarray
) -> tuple[float, float, float]:
    """
    Shoelace area and radial extent of the loop, in the plane basis (e1, e2).

    The loop vertices come out of VTK unordered, so they are sorted by polar
    angle about the centroid before the shoelace sum.  Valid for the
    star-shaped loops a neck orifice produces.
    """
    d = points - origin
    x, y = d @ e1, d @ e2
    order = np.argsort(np.arctan2(y, x))
    xs, ys = x[order], y[order]
    area = 0.5 * abs(np.sum(xs * np.roll(ys, -1) - np.roll(xs, -1) * ys))
    r = np.sqrt(x ** 2 + y ** 2)
    return float(area), float(r.min()), float(r.max())


def extract_neck_orifice(sac_stl: Path) -> NeckPlane:
    """
    Fit the neck plane to the open boundary loop of ``sac_stl``.

    Raises NeckGeometryError if the sac is closed (no open loop), has more than
    one open loop, or the loop is degenerate.
    """
    sac = read_stl(sac_stl)

    edges = sac.extract_feature_edges(
        boundary_edges=True,
        feature_edges=False,
        non_manifold_edges=False,
        manifold_edges=False,
    )
    if edges.n_cells == 0:
        raise NeckGeometryError(
            f"'{Path(sac_stl).name}' is a closed surface — no open boundary loop "
            "found, so the neck orifice cannot be located. VORTEX must clip the "
            "sac at the neck; a capped/watertight sac cannot define an orifice."
        )

    n_regions = int(np.unique(
        edges.connectivity("all").point_data["RegionId"]
    ).size)
    if n_regions != 1:
        raise NeckGeometryError(
            f"found {n_regions} open boundary loops in '{Path(sac_stl).name}'; "
            "expected exactly 1 (the neck orifice). The sac may be clipped at "
            "both ends, or contain holes in the dome."
        )

    pts = np.asarray(edges.points, dtype=float)
    if len(pts) < 3:
        raise NeckGeometryError(
            f"the open boundary loop in '{Path(sac_stl).name}' has only "
            f"{len(pts)} points — too few to define a plane."
        )

    origin = pts.mean(axis=0)
    # Least-variance direction of the loop is the plane normal.
    _, sing, vt = np.linalg.svd(pts - origin)
    normal = vt[2] / np.linalg.norm(vt[2])
    planarity = float(sing[2] / sing[0]) if sing[0] > 0 else 0.0
    if planarity > PLANARITY_WARN:
        warnings.warn(
            f"neck loop in '{Path(sac_stl).name}' is markedly non-planar "
            f"(planarity {planarity:.3g} > {PLANARITY_WARN}); the fitted plane "
            "may not represent the orifice well.",
            stacklevel=2,
        )

    # Orient the normal INTO the sac. Never take the sign from neck_plane.json:
    # a flipped normal would report outflow as inflow, undetectably.
    sac_centroid, centres, areas = _area_weighted_centroid(sac)
    if np.dot(sac_centroid - origin, normal) < 0:
        normal = -normal

    signed = (centres - origin) @ normal
    area_fraction = float(areas[signed > 0].sum() / areas.sum())
    if area_fraction < MIN_AREA_FRACTION:
        warnings.warn(
            f"only {area_fraction:.1%} of the sac area lies on the inward side of "
            f"the fitted neck plane (expected >{MIN_AREA_FRACTION:.0%}). The sac "
            "may not be a dome clipped about the neck, so the normal orientation "
            "is ambiguous.",
            stacklevel=2,
        )

    e1 = vt[0] / np.linalg.norm(vt[0])
    e2 = np.cross(normal, e1)
    loop_area, r_min, r_max = _loop_area_and_radii(pts, origin, e1, e2)
    if not np.isfinite(loop_area) or loop_area <= 0:
        raise NeckGeometryError(
            f"the open boundary loop in '{Path(sac_stl).name}' encloses a "
            f"non-physical area ({loop_area}); the loop is degenerate."
        )
    r_eff = float(np.sqrt(loop_area / np.pi))

    return NeckPlane(
        origin=origin,
        normal=normal,
        radius=r_eff,
        r_min=r_min,
        r_max=r_max,
        r_eff=r_eff,
        loop_area_m2=loop_area,
        planarity=planarity,
        n_loop_points=len(pts),
        area_fraction_positive_side=area_fraction,
    )


def _cross_check(plane: NeckPlane, json_path: Path, origin_scale: float) -> dict:
    """Compare the fitted plane against neck_plane.json (diagnostic only)."""
    data = json.loads(Path(json_path).read_text())
    j_origin = np.asarray(data["origin"], dtype=float) * origin_scale
    j_normal = np.asarray(data["normal"], dtype=float)
    norm = np.linalg.norm(j_normal)
    if norm > 0:
        j_normal = j_normal / norm

    offset = float(np.linalg.norm(j_origin - plane.origin))
    cos = float(np.clip(abs(np.dot(j_normal, plane.normal)), -1.0, 1.0))
    angle = float(np.degrees(np.arccos(cos)))
    return {
        "json_origin_offset_m": offset,
        "json_normal_angle_deg": angle,
        "json_normal_was_unit": bool(abs(norm - 1.0) < 1e-6),
        "json_normal_sign_agrees": bool(np.dot(j_normal, plane.normal) > 0),
        "agrees": bool(offset < max(0.25 * plane.radius, 1e-4) and angle < 15.0),
    }


def resolve_neck_plane(
    sac_stl: Path,
    neck_plane_json: Path | None = None,
    json_origin_scale: float = 1.0,
) -> dict:
    """
    Fit the orifice to ``sac_stl`` and return a JSON-serialisable description.

    ``neck_plane_json`` is compared against the fit and recorded under
    ``cross_check``; it never influences the result.  ``json_origin_scale``
    converts that file's origin into metres (0.001 when the source STLs were mm).
    """
    plane = extract_neck_orifice(sac_stl)
    out = plane.to_dict()
    out["cross_check"] = None
    if neck_plane_json is not None and Path(neck_plane_json).exists():
        try:
            out["cross_check"] = _cross_check(plane, neck_plane_json, json_origin_scale)
        except (KeyError, ValueError, OSError) as e:
            out["cross_check"] = {"error": f"could not read neck_plane.json: {e}"}
    return out


def load_neck_plane(case_dir: Path) -> NeckPlane | None:
    """Read back the resolved plane written into the case directory."""
    path = Path(case_dir) / RESOLVED_FILENAME
    if not path.exists():
        return None
    try:
        data = json.loads(path.read_text())
    except (json.JSONDecodeError, OSError):
        return None
    if "radius_m" not in data:
        return None
    return NeckPlane.from_dict(data)


# ---------------------------------------------------------------------------
# Sampling — slice the volume field at the orifice
# ---------------------------------------------------------------------------

def _cell_vectors(mesh, name: str) -> np.ndarray:
    """Fetch a vector field as cell data, falling back to point data."""
    if name in mesh.cell_data:
        return np.asarray(mesh.cell_data[name], dtype=float)
    if name in mesh.point_data:
        return np.asarray(mesh.point_data[name], dtype=float)
    raise KeyError(
        f"field '{name}' not found on the sliced mesh "
        f"(cell_data={list(mesh.cell_data.keys())}, "
        f"point_data={list(mesh.point_data.keys())})"
    )


def disc_sample(grid, plane: NeckPlane, field: str = "U") -> tuple[np.ndarray, np.ndarray]:
    """
    Slice ``grid`` at the neck plane, clip to the orifice disc, and return
    ``(areas, vectors)`` per cut cell.

    Returns empty arrays (with a warning) when the plane misses the mesh or the
    disc is too small to resolve — never silently reports a zero flux.
    """
    origin = np.asarray(plane.origin, dtype=float)
    normal = np.asarray(plane.normal, dtype=float)

    sl = grid.slice(normal=normal, origin=origin)
    if sl.n_cells == 0:
        warnings.warn(
            "the neck plane does not intersect the mesh — no cells were cut. "
            "Neck metrics are meaningless for this case.",
            stacklevel=2,
        )
        return np.zeros(0), np.zeros((0, 3))

    # In-plane radius: the normal component must be projected out, or oblique
    # planes give the wrong extent.
    d = np.asarray(sl.points, dtype=float) - origin
    radial = d - np.outer(d @ normal, normal)
    sl.point_data["_neck_r"] = np.linalg.norm(radial, axis=1)

    # clip_scalar on point data cuts through cells; a cell-centre mask would
    # give a staircase edge (measured +3.4% area error vs -1.3% for this).
    disc = sl.clip_scalar(scalars="_neck_r", value=float(plane.radius), invert=True)
    if disc.n_cells == 0:
        warnings.warn(
            f"the neck disc (radius {plane.radius:.4g} m) contains no cells — "
            "the orifice is smaller than one mesh cell.",
            stacklevel=2,
        )
        return np.zeros(0), np.zeros((0, 3))
    if disc.n_cells < MIN_DISC_CELLS:
        warnings.warn(
            f"the neck disc contains only {disc.n_cells} cells "
            f"(< {MIN_DISC_CELLS}); the flux quadrature is poorly resolved.",
            stacklevel=2,
        )

    sized = disc.compute_cell_sizes(length=False, area=True, volume=False)
    areas = np.asarray(sized.cell_data["Area"], dtype=float)
    vectors = _cell_vectors(disc, field)
    return areas, vectors


# ---------------------------------------------------------------------------
# Pure metrics — numpy only, no pyvista
# ---------------------------------------------------------------------------

def flux_metrics(areas: np.ndarray, vectors: np.ndarray, normal: np.ndarray) -> dict:
    """
    Flux quantities over the orifice disc for one snapshot.

    inflow_rate  : integral of the INWARD part only, max(U.n, 0) — the clinical
                   "neck inflow rate".  Stays positive even when the net flux
                   through a sealed sac neck averages to ~0.
    net_flux     : signed integral of U.n.  Its cycle mean must be ~0 for a
                   sealed sac; that is the validation instrument for the disc.
    """
    areas = np.asarray(areas, dtype=float)
    vectors = np.asarray(vectors, dtype=float).reshape(-1, 3)
    normal = np.asarray(normal, dtype=float)
    normal = normal / np.linalg.norm(normal)

    if areas.size == 0:
        warnings.warn("empty neck disc — flux metrics are all zero.", stacklevel=2)
        return {
            "inflow_rate_m3s": 0.0,
            "outflow_rate_m3s": 0.0,
            "net_flux_m3s": 0.0,
            "peak_velocity_ms": 0.0,
            "peak_inflow_velocity_ms": 0.0,
            "disc_area_m2": 0.0,
            "n_cells": 0,
        }

    u_n = vectors @ normal
    return {
        "inflow_rate_m3s": float((areas * np.maximum(u_n, 0.0)).sum()),
        "outflow_rate_m3s": float((areas * np.minimum(u_n, 0.0)).sum()),
        "net_flux_m3s": float((areas * u_n).sum()),
        "peak_velocity_ms": float(np.linalg.norm(vectors, axis=1).max()),
        "peak_inflow_velocity_ms": float(max(u_n.max(), 0.0)),
        "disc_area_m2": float(areas.sum()),
        "n_cells": int(areas.size),
    }


# ---------------------------------------------------------------------------
# I/O — walk the solved snapshots
# ---------------------------------------------------------------------------

def read_neck_series(
    case_dir: Path,
    times: list[float],
    plane: NeckPlane,
    field: str = "U",
) -> list[dict]:
    """
    Compute ``flux_metrics`` at each of ``times``.

    Metrics are computed per snapshot and only then time-averaged by the caller;
    they must never be derived from a mean field, since
    ``mean(max(U.n, 0)) != max(mean(U).n, 0)``.
    """
    case_dir = Path(case_dir)
    foam_files = sorted(case_dir.glob("*.foam"))
    if not foam_files:
        raise FileNotFoundError(f"No .foam file in {case_dir}")

    reader = pv.OpenFOAMReader(str(foam_files[0]))
    reader.cell_to_point_creation = False

    rows = []
    for t in times:
        reader.set_active_time_value(float(t))
        mb = reader.read()
        if "internalMesh" not in mb.keys():
            raise KeyError(
                f"'internalMesh' block missing from {foam_files[0].name} at t={t}; "
                "neck metrics need the volume field, not just the boundary patches."
            )
        areas, vectors = disc_sample(mb["internalMesh"], plane, field=field)
        row = flux_metrics(areas, vectors, plane.normal)
        row["time"] = float(t)
        rows.append(row)
    return rows
