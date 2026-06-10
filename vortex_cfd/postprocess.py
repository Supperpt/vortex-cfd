"""
Phase C post-processing: hemodynamic biomarkers from the wall shear stress field.

OpenFOAM's incompressible ``wallShearStress`` function object writes the
**kinematic** WSS vector (tau / rho, units m^2/s^2).  Physical WSS in Pa is
``rho * value``.  OSI is a ratio of WSS integrals and is therefore dimensionless
— the rho factor cancels — so only TAWSS and the Pa thresholds need converting.

Layout:
  * Pure functions (``tawss``, ``osi``, ``summary_stats``) operate on numpy
    arrays only and need no OpenFOAM — they carry the unit tests.
  * The I/O layer (``read_wss_series``, ``compute_metrics``) reads a solved case
    via pyvista's OpenFOAM reader and writes ``metrics_report.json``.

Authoritative TAWSS / OSI come from here, computed from the per-snapshot
wallShearStress fields over the last cardiac cycle.  The ``fieldAverage``
function object's ``wallShearStressMean`` is only a ParaView convenience.

New (default) mode: metrics are scoped to the ``aneurysm_sac`` patch.
  Normalised WSS = mean sac TAWSS / mean parent_vessel TAWSS.
  Neck inflow rate and sac pressure are parsed from postProcessing/ CSVs.

Legacy mode: metrics are scoped to the single ``wall`` patch as in Phase C.
  Mode is auto-detected from ``patch_labels.json`` in the case directory.
"""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pyvista as pv

from .waveform import T_CYCLE

# Default fluid properties (kept in sync with case_builder defaults).
RHO = 1060.0      # kg/m^3, whole blood
NU = 3.3e-6       # m^2/s, kinematic viscosity

# Validation ranges (warn — never abort — if a metric falls outside).
WSS_MAX_PA = 50.0
OSI_MAX = 0.5
# Clinical "low-and-oscillatory" risk thresholds.
TAWSS_LOW_PA = 0.4
OSI_HIGH = 0.3


# ---------------------------------------------------------------------------
# Pure functions — numpy only, no OpenFOAM
# ---------------------------------------------------------------------------

def _time_weights(times: np.ndarray) -> np.ndarray:
    """
    Trapezoidal integration weights for samples taken at ``times``.
    The weights sum to (times[-1] - times[0]); a uniformly-spaced series
    reduces to equal weights, i.e. a plain time average.
    """
    times = np.asarray(times, dtype=float)
    if times.size == 1:
        return np.array([1.0])
    w = np.zeros_like(times)
    w[1:-1] = 0.5 * (times[2:] - times[:-2])
    w[0] = 0.5 * (times[1] - times[0])
    w[-1] = 0.5 * (times[-1] - times[-2])
    return w


def tawss(wss_series: np.ndarray, weights: np.ndarray | None = None) -> np.ndarray:
    """
    Time-averaged WSS magnitude (kinematic) per wall face.

    wss_series : (n_times, n_faces, 3) kinematic WSS vectors.
    weights    : optional (n_times,) integration weights; uniform if None.
    Returns    : (n_faces,) array.
    """
    mag = np.linalg.norm(wss_series, axis=2)            # (n_times, n_faces)
    if weights is None:
        return mag.mean(axis=0)
    weights = np.asarray(weights, dtype=float)
    return (weights[:, None] * mag).sum(axis=0) / weights.sum()


def osi(wss_series: np.ndarray, weights: np.ndarray | None = None) -> np.ndarray:
    """
    Oscillatory Shear Index per wall face:

        OSI = 0.5 * (1 - |mean(WSS_vec)| / mean(|WSS|))

    Dimensionless (0 = unidirectional flow, 0.5 = fully reversing).
    Faces with zero mean magnitude are reported as OSI = 0.
    """
    if weights is None:
        mean_vec = wss_series.mean(axis=0)                       # (n_faces, 3)
        mean_mag = np.linalg.norm(wss_series, axis=2).mean(axis=0)  # (n_faces,)
    else:
        weights = np.asarray(weights, dtype=float)
        wsum = weights.sum()
        mean_vec = (weights[:, None, None] * wss_series).sum(axis=0) / wsum
        mean_mag = (weights[:, None] * np.linalg.norm(wss_series, axis=2)).sum(axis=0) / wsum

    vec_mag = np.linalg.norm(mean_vec, axis=1)                   # (n_faces,)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(mean_mag > 0, vec_mag / mean_mag, 1.0)
    return 0.5 * (1.0 - ratio)


def summary_stats(
    tawss_kin: np.ndarray,
    osi_field: np.ndarray,
    areas: np.ndarray,
    rho: float = RHO,
) -> dict:
    """
    Area-weighted summary statistics over the wall.  All means are weighted by
    face area (faces vary in size); area fractions use the same areas.
    """
    areas = np.asarray(areas, dtype=float)
    total_area = float(areas.sum())
    tawss_pa = rho * tawss_kin

    def area_fraction(mask: np.ndarray) -> float:
        if total_area <= 0:
            return 0.0
        return float(areas[mask].sum() / total_area)

    def wmean(field: np.ndarray) -> float:
        if total_area <= 0:
            return float(field.mean())
        return float(np.average(field, weights=areas))

    return {
        "tawss_pa": {
            "mean": wmean(tawss_pa),
            "max": float(tawss_pa.max()),
            "min": float(tawss_pa.min()),
        },
        "tawss_kinematic": {
            "mean": wmean(tawss_kin),
            "max": float(tawss_kin.max()),
        },
        "osi": {
            "mean": wmean(osi_field),
            "max": float(osi_field.max()),
        },
        "area_fraction_osi_gt_0p3": area_fraction(osi_field > OSI_HIGH),
        "area_fraction_tawss_lt_0p4pa": area_fraction(tawss_pa < TAWSS_LOW_PA),
        "wall_area_m2": total_area,
        "n_wall_faces": int(tawss_kin.size),
    }


# ---------------------------------------------------------------------------
# I/O layer — reads a solved OpenFOAM case via pyvista
# ---------------------------------------------------------------------------

def _wall_block(multiblock, wall_patch: str):
    """Locate the wall boundary patch inside the OpenFOAM reader output."""
    # ESI cases expose patches under a 'boundary' MultiBlock.
    if "boundary" in multiblock.keys():
        boundary = multiblock["boundary"]
        if wall_patch in boundary.keys():
            return boundary[wall_patch]
    # Fallback: search every nested block for one carrying wallShearStress.
    for name in multiblock.keys():
        block = multiblock[name]
        if block is None:
            continue
        if hasattr(block, "keys"):
            found = _wall_block(block, wall_patch)
            if found is not None:
                return found
        elif "wallShearStress" in getattr(block, "cell_data", {}):
            return block
    return None


def _wss_from_block(block) -> np.ndarray:
    """Return the (n_faces, 3) wallShearStress array (cell data preferred)."""
    if "wallShearStress" in block.cell_data:
        return np.asarray(block.cell_data["wallShearStress"])
    if "wallShearStress" in block.point_data:
        return np.asarray(block.point_data["wallShearStress"])
    raise KeyError("wallShearStress not present on the wall patch — was the "
                   "function object enabled (--postprocess) or generated via "
                   "-postProcess?")


def available_times(case_dir: Path) -> list[float]:
    """Snapshot times present in the case, read from .foam metadata only."""
    case_dir = Path(case_dir)
    foam_files = sorted(case_dir.glob("*.foam"))
    if not foam_files:
        raise FileNotFoundError(f"No .foam file in {case_dir}")
    reader = pv.OpenFOAMReader(str(foam_files[0]))
    return [float(t) for t in reader.time_values]


def read_wss_series(
    case_dir: Path,
    t_start: float,
    wall_patch: str = "wall",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Read the wall wallShearStress field at every snapshot with time >= t_start.

    Returns (times, wss, areas):
        times : (n_times,)            absolute snapshot times [s]
        wss   : (n_times, n_faces, 3) kinematic WSS vectors
        areas : (n_faces,)            wall face areas [m^2]
    """
    case_dir = Path(case_dir)
    foam_files = sorted(case_dir.glob("*.foam"))
    if not foam_files:
        raise FileNotFoundError(f"No .foam file in {case_dir}")

    reader = pv.OpenFOAMReader(str(foam_files[0]))
    reader.enable_all_patch_arrays()
    # Keep face (cell) values rather than interpolating to points.
    reader.cell_to_point_creation = False

    times = [t for t in reader.time_values if t >= t_start - 1e-9]
    if not times:
        raise ValueError(
            f"No snapshots at or after t_start={t_start:.4f}s "
            f"(available: {list(reader.time_values)})"
        )

    series = []
    areas = None
    for t in times:
        reader.set_active_time_value(t)
        mb = reader.read()
        wall = _wall_block(mb, wall_patch)
        if wall is None:
            raise KeyError(f"Wall patch '{wall_patch}' not found in {foam_files[0]}")
        series.append(_wss_from_block(wall))
        if areas is None:
            sized = wall.compute_cell_sizes(length=False, area=True, volume=False)
            areas = np.asarray(sized.cell_data["Area"])

    return np.asarray(times, dtype=float), np.asarray(series), areas


def _detect_patches(case_dir: Path) -> tuple[str | None, str]:
    """
    Read patch_labels.json to detect pipeline mode.
    Returns (aneurysm_patch, parent_vessel_patch):
      - New mode: ("aneurysm_sac", "parent_vessel")
      - Legacy mode: (None, "wall")
    """
    labels_path = Path(case_dir) / "patch_labels.json"
    if labels_path.exists():
        labels = json.loads(labels_path.read_text())
        if "aneurysm_sac" in labels:
            return "aneurysm_sac", "parent_vessel"
    return None, "wall"


def _read_surface_field_value(
    case_dir: Path,
    fo_name: str,
    t_start: float,
) -> list[tuple[float, float]]:
    """
    Parse a surfaceFieldValue postProcessing output file.

    Reads from postProcessing/<fo_name>/<startTime>/surfaceFieldValue.dat.
    Returns [(time, value), ...] for times >= t_start.

    For scalar fields: value is the scalar.
    For vector fields (3 components after time column): value is the magnitude.
    Lines beginning with '#' are skipped as comments.
    """
    pp_dir = Path(case_dir) / "postProcessing" / fo_name
    if not pp_dir.exists():
        return []

    # Use the earliest time directory (matches the simulation startTime).
    time_dirs = sorted(pp_dir.iterdir())
    if not time_dirs:
        return []

    dat_path = time_dirs[0] / "surfaceFieldValue.dat"
    if not dat_path.exists():
        return []

    rows: list[tuple[float, float]] = []
    with dat_path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                t = float(parts[0])
                if t < t_start - 1e-9:
                    continue
                if len(parts) == 2:
                    rows.append((t, float(parts[1])))
                elif len(parts) == 4:
                    # Vector: (time, vx, vy, vz) — report magnitude.
                    vx, vy, vz = float(parts[1]), float(parts[2]), float(parts[3])
                    rows.append((t, float(np.sqrt(vx**2 + vy**2 + vz**2))))
            except ValueError:
                pass
    return rows


# ---------------------------------------------------------------------------
# Orchestration — read + math + report
# ---------------------------------------------------------------------------

def compute_metrics(
    case_dir: Path,
    cycles: int | None = None,
    t_cycle: float = T_CYCLE,
    rho: float = RHO,
    nu: float = NU,
    wall_patch: str = "wall",
) -> dict:
    """
    Compute TAWSS / OSI over the last cardiac cycle and write
    ``metrics_report.json`` into the case directory.  Returns the report dict.

    Mode is auto-detected from ``patch_labels.json``:
      - New mode (aneurysm_sac present): metrics scoped to aneurysm_sac patch;
        also computes parent_vessel mean TAWSS and normalised WSS.
      - Legacy mode: metrics scoped to the wall patch (Phase C behaviour).

    If ``cycles`` is given the last cycle is [(cycles-1)*T, cycles*T].  If it is
    None (standalone post-processing of an old case where the cycle count is
    unknown) the last cycle is derived from the data: the final available
    snapshot time minus one period.

    Validation ranges are checked and warned about — never aborts.
    """
    case_dir = Path(case_dir)
    if cycles is not None:
        t_start = (cycles - 1) * t_cycle
    else:
        t_max = max(available_times(case_dir))
        t_start = max(0.0, t_max - t_cycle)

    aneurysm_patch, parent_vessel_patch = _detect_patches(case_dir)
    is_new_mode = aneurysm_patch is not None

    if is_new_mode:
        print(f"\n[vortex-cfd] Post-processing WSS/TAWSS/OSI (aneurysm sac) "
              f"over the last cycle (t >= {t_start:.4f}s) ...")
    else:
        print(f"\n[vortex-cfd] Post-processing WSS/TAWSS/OSI over the last cycle "
              f"(t >= {t_start:.4f}s) ...")

    # Primary patch: aneurysm_sac in new mode, wall in legacy mode.
    primary_patch = aneurysm_patch if is_new_mode else wall_patch
    times, wss, areas = read_wss_series(case_dir, t_start, primary_patch)
    weights = _time_weights(times)

    tawss_kin = tawss(wss, weights)
    osi_field = osi(wss, weights)
    stats = summary_stats(tawss_kin, osi_field, areas, rho)

    in_wss_range = bool(stats["tawss_pa"]["max"] <= WSS_MAX_PA
                        and stats["tawss_pa"]["min"] >= 0.0)
    in_osi_range = bool(0.0 <= stats["osi"]["max"] <= OSI_MAX)

    report: dict = {
        **stats,
        "cycle_analysed": [float(times[0]), float(times[-1])],
        "n_snapshots": int(times.size),
        "rho": rho,
        "nu": nu,
        "validation": {
            "tawss_pa_in_0_50": in_wss_range,
            "osi_in_0_0.5": in_osi_range,
        },
    }

    if is_new_mode:
        # Rename generic keys to aneurysm-scoped names for clarity.
        report["aneurysm_tawss_pa"] = report.pop("tawss_pa")
        report["aneurysm_tawss_kinematic"] = report.pop("tawss_kinematic")
        report["aneurysm_osi"] = report.pop("osi")
        report["aneurysm_wall_area_m2"] = report.pop("wall_area_m2")
        report["aneurysm_n_faces"] = report.pop("n_wall_faces")

        # Parent vessel mean TAWSS (used as normalisation denominator).
        _, wss_pv, areas_pv = read_wss_series(case_dir, t_start, parent_vessel_patch)
        w_pv = _time_weights(times)
        tawss_pv_kin = tawss(wss_pv, w_pv)
        parent_tawss_pa = float(rho * np.average(tawss_pv_kin,
                                                   weights=np.asarray(areas_pv)))
        report["parent_tawss_pa_mean"] = parent_tawss_pa

        sac_tawss_mean = report["aneurysm_tawss_pa"]["mean"]
        if parent_tawss_pa > 0:
            report["normalised_wss"] = float(sac_tawss_mean / parent_tawss_pa)
        else:
            report["normalised_wss"] = None

        # Sac pressure from surfaceFieldValue postProcessing CSVs.
        press_mean_rows = _read_surface_field_value(
            case_dir, "surfaceFieldValue_sac_pressure_mean", t_start)
        press_max_rows = _read_surface_field_value(
            case_dir, "surfaceFieldValue_sac_pressure_max", t_start)

        if press_mean_rows:
            # Kinematic pressure → Pa
            report["sac_pressure_mean_pa"] = float(
                rho * np.mean([v for _, v in press_mean_rows]))
        else:
            report["sac_pressure_mean_pa"] = None

        if press_max_rows:
            report["sac_pressure_peak_pa"] = float(
                rho * max(v for _, v in press_max_rows))
        else:
            report["sac_pressure_peak_pa"] = None

        # Neck inflow rate from surfaceFieldValue postProcessing CSVs.
        flux_rows = _read_surface_field_value(
            case_dir, "surfaceFieldValue_neck_flux", t_start)
        peak_vel_rows = _read_surface_field_value(
            case_dir, "surfaceFieldValue_neck_peak_vel", t_start)

        if flux_rows:
            vals = [v for _, v in flux_rows]
            report["neck_mean_flow_rate_m3s"] = float(np.mean(vals))
            report["neck_peak_flow_rate_m3s"] = float(max(vals))
        else:
            report["neck_mean_flow_rate_m3s"] = None
            report["neck_peak_flow_rate_m3s"] = None

        if peak_vel_rows:
            report["neck_peak_velocity_ms"] = float(max(v for _, v in peak_vel_rows))
        else:
            report["neck_peak_velocity_ms"] = None

    report["generated"] = datetime.now().isoformat(timespec="seconds")

    out = case_dir / "metrics_report.json"
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")

    # Human-readable summary.
    if is_new_mode:
        at = report["aneurysm_tawss_pa"]
        ao = report["aneurysm_osi"]
        print(f"  Aneurysm faces      : {report['aneurysm_n_faces']} "
              f"({report['aneurysm_wall_area_m2']*1e6:.1f} mm^2)")
        print(f"  Snapshots           : {report['n_snapshots']} "
              f"over [{times[0]:.4f}, {times[-1]:.4f}] s")
        print(f"  TAWSS sac (Pa)      : mean {at['mean']:.3f}, max {at['max']:.3f}")
        print(f"  OSI sac             : mean {ao['mean']:.4f}, max {ao['max']:.4f}")
        print(f"  Parent TAWSS (Pa)   : mean {parent_tawss_pa:.3f}")
        if report["normalised_wss"] is not None:
            print(f"  Normalised WSS      : {report['normalised_wss']:.3f}")
        print(f"  Low-WSS area (<0.4 Pa)   : "
              f"{report['area_fraction_tawss_lt_0p4pa']*100:.1f} %")
        print(f"  High-OSI area (>0.3)     : "
              f"{report['area_fraction_osi_gt_0p3']*100:.1f} %")
    else:
        print(f"  Wall faces analysed : {stats['n_wall_faces']} "
              f"({stats['wall_area_m2']*1e6:.1f} mm^2)")
        print(f"  Snapshots           : {report['n_snapshots']} "
              f"over [{times[0]:.4f}, {times[-1]:.4f}] s")
        print(f"  TAWSS (Pa)          : mean {stats['tawss_pa']['mean']:.3f}, "
              f"max {stats['tawss_pa']['max']:.3f}")
        print(f"  OSI                 : mean {stats['osi']['mean']:.4f}, "
              f"max {stats['osi']['max']:.4f}")
        print(f"  Low-WSS area (<0.4 Pa)   : "
              f"{stats['area_fraction_tawss_lt_0p4pa']*100:.1f} %")
        print(f"  High-OSI area (>0.3)     : "
              f"{stats['area_fraction_osi_gt_0p3']*100:.1f} %")

    if not in_wss_range:
        print(f"  WARNING: TAWSS outside the expected 0–{WSS_MAX_PA:.0f} Pa range.")
    if not in_osi_range:
        print(f"  WARNING: OSI outside the expected 0–{OSI_MAX} range.")

    print(f"  Report written      : {out}")
    return report
