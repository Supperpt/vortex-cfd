"""Build the OpenFOAM case directory from scaled STLs and simulation parameters."""
from __future__ import annotations

import json
import multiprocessing
import shutil
from datetime import datetime
from pathlib import Path

import numpy as np
from jinja2 import Environment, FileSystemLoader

from . import neck
from .waveform import T_CYCLE
from .scaling import _any_in_mm, read_stl
from .constants import RHO, NU

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

def _bbox_with_buffer(wall_stls: "Path | list[Path]", buffer: float = 0.20) -> dict:
    """
    Bounding box of one or more wall STL(s), expanded by `buffer` on each side.
    Accepts a single Path or a list; a list unions the bboxes of all surfaces.
    """
    if isinstance(wall_stls, Path):
        wall_stls = [wall_stls]
    bounds = [read_stl(p).bounds for p in wall_stls]
    xmin = min(b[0] for b in bounds)
    xmax = max(b[1] for b in bounds)
    ymin = min(b[2] for b in bounds)
    ymax = max(b[3] for b in bounds)
    zmin = min(b[4] for b in bounds)
    zmax = max(b[5] for b in bounds)
    dx, dy, dz = xmax - xmin, ymax - ymin, zmax - zmin
    return {
        "xmin": xmin - buffer * dx, "xmax": xmax + buffer * dx,
        "ymin": ymin - buffer * dy, "ymax": ymax + buffer * dy,
        "zmin": zmin - buffer * dz, "zmax": zmax + buffer * dz,
    }


def _background_cell_counts(bbox: dict, target: float = 0.002) -> dict:
    def n(lo, hi):
        return max(10, min(60, round((hi - lo) / target)))
    return {
        "nx": n(bbox["xmin"], bbox["xmax"]),
        "ny": n(bbox["ymin"], bbox["ymax"]),
        "nz": n(bbox["zmin"], bbox["zmax"]),
    }


def _location_in_mesh(inlet_stl: Path) -> tuple[float, float, float]:
    """
    Interior point guaranteed to be inside the lumen: inlet face centroid
    displaced one inlet-radius inward along the area-weighted face normal.

    This replaces the bounding-box centre of the wall STL, which fails for
    curved vessels where the bbox centre falls inside the wall material.
    """
    mesh = read_stl(inlet_stl)
    sized = mesh.compute_cell_sizes()
    areas = sized.cell_data["Area"]
    total_area = areas.sum()
    radius = float(np.sqrt(total_area / np.pi))

    # Area-weighted centroid and normal
    pts = np.array(mesh.cell_centers().points)
    centroid = (areas[:, None] * pts).sum(axis=0) / total_area

    normals = np.array(mesh.compute_normals(cell_normals=True, point_normals=False).cell_data["Normals"])
    avg_normal = (areas[:, None] * normals).sum(axis=0)
    avg_normal /= np.linalg.norm(avg_normal)

    # Step one radius inward (negate normal — cap normals point outward)
    interior = centroid - avg_normal * radius

    return (float(interior[0]), float(interior[1]), float(interior[2]))


def _inlet_area(inlet_stl: Path) -> float:
    mesh = read_stl(inlet_stl)
    sized = mesh.compute_cell_sizes()
    area = float(sized.cell_data["Area"].sum())
    if not np.isfinite(area) or area <= 0:
        raise ValueError(
            f"Inlet STL '{inlet_stl}' produced a non-physical area {area}. "
            "The surface may be corrupt, empty, or degenerate."
        )
    return area


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
    legacy: bool = False,
    stl_source_dir: Path | None = None,
) -> Path:
    """
    Render all Jinja2 templates and assemble the OpenFOAM case directory.
    Returns the Path to the created directory.

    New mode (default): expects aneurysm_sac + parent_vessel in scaled_stls.
    Legacy mode (legacy=True or wall key present): uses a single wall patch.

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
    inlet_stl = case_dir / "constant" / "triSurface" / "inlet.stl"

    # Detect mode: legacy uses a single "wall" patch; new mode uses two named patches.
    is_legacy = legacy or ("wall" in scaled_stls)

    if is_legacy:
        wall_stl = case_dir / "constant" / "triSurface" / "wall.stl"
        bbox = _bbox_with_buffer(wall_stl)
        wall_patches = ["wall"]
        aneurysm_patch = None
        parent_vessel_patch = "wall"
        has_neck_plane = False
        neck_origin: list[float] = [0.0, 0.0, 0.0]
        neck_normal: list[float] = [0.0, 0.0, 1.0]
        neck_resolved = None  # legacy mode has no aneurysm sac
        patch_labels = {"wall": "wall", "inlet": "inlet"}
    else:
        sac_stl = case_dir / "constant" / "triSurface" / "aneurysm_sac.stl"
        pv_stl = case_dir / "constant" / "triSurface" / "parent_vessel.stl"
        bbox = _bbox_with_buffer([sac_stl, pv_stl])
        wall_patches = ["aneurysm_sac", "parent_vessel"]
        aneurysm_patch = "aneurysm_sac"
        parent_vessel_patch = "parent_vessel"
        json_origin_scale = 1.0

        # Load neck_plane.json (or output_neck_plane.json) — skip if absent.
        neck_plane_path = None
        if stl_source_dir:
            for candidate in ("neck_plane.json", "output_neck_plane.json"):
                p = stl_source_dir / candidate
                if p.exists():
                    neck_plane_path = p
                    break
        if neck_plane_path:
            data = json.loads(neck_plane_path.read_text())
            neck_origin = list(data["origin"])
            neck_normal = list(data["normal"])
            # Scale origin mm→m if the source STLs are in millimetres.
            src_stls = list(stl_source_dir.glob("*.stl"))
            if src_stls and _any_in_mm(src_stls):
                json_origin_scale = 0.001
                neck_origin = [v * 0.001 for v in neck_origin]
            has_neck_plane = True
        else:
            if stl_source_dir:
                print(f"WARNING: neck_plane.json not found in {stl_source_dir}. "
                      "Neck-plane function objects will be skipped.")
            neck_origin = [0.0, 0.0, 0.0]
            neck_normal = [0.0, 0.0, 1.0]
            has_neck_plane = False

        # Fit the neck orifice to the sac STL's open boundary loop.  The copy in
        # the case dir is already in metres, so no unit heuristic is needed here.
        # This is a diagnostic: a sac that cannot be fitted must not fail a build.
        try:
            neck_resolved = neck.resolve_neck_plane(
                sac_stl, neck_plane_path, json_origin_scale
            )
        except neck.NeckGeometryError as e:
            print(f"WARNING: neck orifice could not be determined — {e}")
            neck_resolved = {"status": "unavailable", "reason": str(e)}

        patch_labels = {"aneurysm_sac": "aneurysm_sac", "parent_vessel": "parent_vessel",
                        "inlet": "inlet"}

    patch_labels.update({n: "outlet" for n in outlet_names})

    cell_counts = _background_cell_counts(bbox)
    loc = _location_in_mesh(inlet_stl)
    area = _inlet_area(inlet_stl)
    table = _waveform_table(waveform, mean_velocity, area)

    end_time = cycles * T_CYCLE
    write_interval = T_CYCLE / 50  # 50 snapshots per cycle
    # fieldAverage / Python post-processing analyse only the last cycle; earlier
    # cycles are discarded as transient initialisation.
    field_average_start = (cycles - 1) * T_CYCLE

    ctx = {
        "wall_patches":        wall_patches,
        "aneurysm_patch":      aneurysm_patch,
        "parent_vessel_patch": parent_vessel_patch,
        "has_neck_plane":      has_neck_plane,
        "neck_origin":         neck_origin,
        "neck_normal":         neck_normal,
        "inlet_patch":         "inlet",
        "outlet_patches":      outlet_names,
        "all_stls":            list(scaled_stls.keys()),
        "waveform_table":      table,
        "t_cycle":             T_CYCLE,
        "end_time":            end_time,
        "write_interval":      write_interval,
        "max_co":              0.8,
        "cores":               cores,
        "nu":                  NU,
        "rho":                 RHO,
        "bbox":                bbox,
        "nx":                  cell_counts["nx"],
        "ny":                  cell_counts["ny"],
        "nz":                  cell_counts["nz"],
        "location_in_mesh":    loc,
        "case_name":           case_name,
        "postprocess":         postprocess,
        "field_average_start": field_average_start,
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
    (case_dir / "patch_labels.json").write_text(
        json.dumps(patch_labels, indent=2), encoding="utf-8"
    )

    # The fitted neck orifice, so post-processing need not re-read the STLs.
    if neck_resolved is not None:
        (case_dir / neck.RESOLVED_FILENAME).write_text(
            json.dumps(neck_resolved, indent=2), encoding="utf-8"
        )

    # ParaView placeholder
    (case_dir / f"{case_name}.foam").touch()

    return case_dir
