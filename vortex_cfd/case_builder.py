"""Build the OpenFOAM case directory from scaled STLs and simulation parameters."""

import json
import multiprocessing
import shutil
from datetime import datetime
from pathlib import Path

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


def _location_in_mesh(wall_stl: Path) -> tuple[float, float, float]:
    """Approximate interior point: centre of the wall STL bounding box."""
    mesh = pv.read(str(wall_stl))
    c = mesh.center
    return (float(c[0]), float(c[1]), float(c[2]))


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
) -> Path:
    """
    Render all Jinja2 templates and assemble the OpenFOAM case directory.
    Returns the Path to the created directory.

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
    loc = _location_in_mesh(wall_stl)
    area = _inlet_area(inlet_stl)
    table = _waveform_table(waveform, mean_velocity, area)

    end_time = cycles * T_CYCLE
    write_interval = T_CYCLE / 50  # 50 snapshots per cycle

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

    return case_dir
