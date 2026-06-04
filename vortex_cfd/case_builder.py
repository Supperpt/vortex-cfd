"""Build the OpenFOAM case directory from scaled STLs and simulation parameters."""

import json
import logging
import multiprocessing
import shutil
from datetime import datetime
from pathlib import Path

log = logging.getLogger("vortex_cfd")

import numpy as np
import pyvista as pv
from jinja2 import Environment, FileSystemLoader

from .waveform import T_CYCLE

TEMPLATES_DIR = Path(__file__).parent / "templates"

_TEMPLATE_MAP = {
    "0/U.j2":                              "0/U",
    "0/p.j2":                              "0/p",
    "constant/transportProperties.j2":    "constant/transportProperties",
    "constant/turbulenceProperties.j2":   "constant/turbulenceProperties",
    "system/controlDict.j2":              "system/controlDict",
    "system/fvSchemes.j2":                "system/fvSchemes",
    "system/fvSolution.j2":               "system/fvSolution",
    "system/snappyHexMeshDict.j2":        "system/snappyHexMeshDict",
    "system/decomposeParDict.j2":         "system/decomposeParDict",
    "system/meshQualityDict.j2":          "system/meshQualityDict",
    "system/surfaceFeatureExtractDict.j2": "system/surfaceFeatureExtractDict",
    "system/blockMeshDict.j2":            "system/blockMeshDict",
    "Allrun.j2":                          "Allrun",
}


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def _bbox_with_buffer(wall_stl: Path, buffer: float = 0.20) -> dict:
    mesh = pv.read(str(wall_stl))
    b = mesh.bounds  # (xmin, xmax, ymin, ymax, zmin, zmax)
    dx, dy, dz = b[1] - b[0], b[3] - b[2], b[5] - b[4]
    return {
        "xmin": b[0] - buffer * dx, "xmax": b[1] + buffer * dx,
        "ymin": b[2] - buffer * dy, "ymax": b[3] + buffer * dy,
        "zmin": b[4] - buffer * dz, "zmax": b[5] + buffer * dz,
    }


def _background_cell_counts(bbox: dict, target: float = 0.002) -> dict:
    def n(lo, hi):
        return max(10, min(60, round((hi - lo) / target)))
    return {
        "nx": n(bbox["xmin"], bbox["xmax"]),
        "ny": n(bbox["ymin"], bbox["ymax"]),
        "nz": n(bbox["zmin"], bbox["zmax"]),
    }


def _inlet_geometry(inlet_stl: Path) -> tuple:
    """
    Compute the geometry of the inlet cap STL:
        centroid : (3,) area-weighted face centroid [m]
        normal   : (3,) unit inward normal (VMTK cap normals point outward; negated here)
        radius   : float, equivalent circular radius [m]
        interior : (3,) point one radius inward — used as locationInMesh

    The interior point is guaranteed inside the lumen for curved vessels.
    See BUG-006 and D-003.  Replaces the original bounding-box centre approach.
    """
    mesh = pv.read(str(inlet_stl))
    sized = mesh.compute_cell_sizes()
    areas = sized.cell_data["Area"]
    total_area = areas.sum()
    radius = float(np.sqrt(total_area / np.pi))

    pts = np.array(mesh.cell_centers().points)
    centroid = (areas[:, None] * pts).sum(axis=0) / total_area

    normals = np.array(mesh.compute_normals(cell_normals=True, point_normals=False).cell_data["Normals"])
    avg_normal = (areas[:, None] * normals).sum(axis=0)
    avg_normal /= np.linalg.norm(avg_normal)

    inward_normal = -avg_normal          # cap normals point outward by VMTK convention
    interior = centroid + inward_normal * radius
    return (
        centroid.astype(float),
        inward_normal.astype(float),
        radius,
        (float(interior[0]), float(interior[1]), float(interior[2])),
    )


def _location_in_mesh(inlet_stl: Path) -> tuple[float, float, float]:
    """Thin wrapper — backward-compatible entry point used by existing tests."""
    _, _, _, interior = _inlet_geometry(inlet_stl)
    return interior


def _inlet_area(inlet_stl: Path) -> float:
    mesh = pv.read(str(inlet_stl))
    sized = mesh.compute_cell_sizes()
    return float(sized.cell_data["Area"].sum())


def _waveform_table(
    waveform: np.ndarray,
    mean_velocity: float,
    inlet_area: float,
    t_cycle: float = T_CYCLE,
) -> list[tuple[float, float]]:
    Q_mean = mean_velocity * inlet_area
    t_abs = waveform[:, 0] * t_cycle
    Q_abs = waveform[:, 1] * Q_mean
    return list(zip(t_abs.tolist(), Q_abs.tolist()))


# ---------------------------------------------------------------------------
# Case builder
# ---------------------------------------------------------------------------

def build_case(
    scaled_stls: dict[str, Path],
    labels: dict,
    cycles: int,
    mean_velocity: float,
    waveform: np.ndarray,
    cores: int | None,
    out_dir: str,
    postprocess: bool = False,
    womersley: bool = False,
) -> tuple:
    """
    Render all Jinja2 templates and assemble the OpenFOAM case directory.

    Returns (case_dir, inlet_params) where inlet_params is a dict with keys
    'centroid', 'normal', 'radius', 'mean_velocity', 'waveform', 'nu'
    needed by runner.py to generate Womersley boundary data after meshing.

    patch_labels.json is written after the directory exists, not before — avoiding
    the race condition that plagued the prior VesselForge_AutoCFD iteration.
    """
    if cores is None:
        cores = multiprocessing.cpu_count()

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    case_name = f"case_{timestamp}"
    case_dir = Path(out_dir).resolve() / case_name

    # Build directory tree first
    for sub in ["0", "constant/triSurface", "system"]:
        (case_dir / sub).mkdir(parents=True)

    # Copy scaled STLs — must happen before patch_labels.json is written
    for canonical, src in scaled_stls.items():
        shutil.copy2(str(src), str(case_dir / "constant" / "triSurface" / f"{canonical}.stl"))

    outlet_names = sorted(k for k in scaled_stls if k.startswith("outlet_"))
    wall_stl = case_dir / "constant" / "triSurface" / "wall.stl"
    inlet_stl = case_dir / "constant" / "triSurface" / "inlet.stl"

    bbox = _bbox_with_buffer(wall_stl)
    cell_counts = _background_cell_counts(bbox)
    centroid, normal, radius, loc = _inlet_geometry(inlet_stl)
    area = float(np.pi * radius ** 2)
    table = _waveform_table(waveform, mean_velocity, area)

    # Physiological flow-rate validation (warn only — never abort)
    Q_mL = mean_velocity * area * 1e6
    if not (1.0 <= Q_mL <= 10.0):
        log.warning("Mean flow rate %.2f mL/s is outside the typical ICA range "
                    "(1–10 mL/s). Check --mean-velocity and ensure the STL is in metres.",
                    Q_mL)

    end_time = cycles * T_CYCLE
    write_interval = T_CYCLE / 50  # 50 snapshots per cycle
    # fieldAverage / Python post-processing analyse only the last cycle; earlier
    # cycles are discarded as transient initialisation.
    field_average_start = (cycles - 1) * T_CYCLE

    ctx = {
        "wall_patch":      "wall",
        "inlet_patch":     "inlet",
        "outlet_patches":  outlet_names,
        "all_stls":        list(scaled_stls.keys()),
        "waveform_table":  table,
        "t_cycle":         T_CYCLE,
        "end_time":        end_time,
        "write_interval":  write_interval,
        "max_co":          0.8,
        "cores":           cores,
        "nu":              3.3e-6,
        "rho":             1060.0,
        "bbox":            bbox,
        "nx":              cell_counts["nx"],
        "ny":              cell_counts["ny"],
        "nz":              cell_counts["nz"],
        "location_in_mesh": loc,
        "case_name":       case_name,
        "postprocess":     postprocess,
        "field_average_start": field_average_start,
        "womersley":       womersley,
    }

    jinja_env = Environment(
        loader=FileSystemLoader(str(TEMPLATES_DIR)),
        keep_trailing_newline=True,
    )

    for tmpl_name, dest_rel in _TEMPLATE_MAP.items():
        rendered = jinja_env.get_template(tmpl_name).render(**ctx)
        dest = case_dir / dest_rel
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(rendered, encoding="utf-8")

    # chmod +x Allrun
    allrun = case_dir / "Allrun"
    allrun.chmod(allrun.stat().st_mode | 0o111)

    # Metadata — written after the directory exists (lesson from prior iteration)
    patch_labels = {"wall": "wall", "inlet": "inlet"}
    patch_labels.update({n: "outlet" for n in outlet_names})
    (case_dir / "patch_labels.json").write_text(
        json.dumps(patch_labels, indent=2), encoding="utf-8"
    )

    # ParaView placeholder
    (case_dir / f"{case_name}.foam").touch()

    # Params forwarded to runner.py:
    #   inlet_params  — Womersley geometry + flow data for post-mesh boundaryData
    #   mesh_params   — bbox + base Jinja context for SHM retry blockMesh re-render
    inlet_params = {
        "centroid":      centroid,
        "normal":        normal,
        "radius":        radius,
        "mean_velocity": mean_velocity,
        "waveform":      waveform,
        "nu":            3.3e-6,
        # SHM retry needs to re-render blockMeshDict with a finer target
        "_bbox":         bbox,
        "_ctx":          ctx,
        "_tmpl_dir":     str(TEMPLATES_DIR),
    }

    return case_dir, inlet_params
