"""Orchestrate the full OpenFOAM meshing and solving pipeline."""

import logging
import multiprocessing
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

from . import postprocess

log = logging.getLogger("vortex_cfd")


def _build_env(of_env: dict) -> dict | None:
    """Merge captured OpenFOAM env vars into the current process environment."""
    extra = of_env.get("env")
    if not extra:
        return None  # OpenFOAM already sourced; inherit as-is
    merged = dict(os.environ)
    merged.update(extra)
    return merged


def _run(cmd: str, cwd: Path, env: dict | None, label: str) -> str:
    """Run cmd, stream output to stdout, return combined stdout+stderr as string."""
    log.info("\n[vortex-cfd] %s", label)
    log.info("  $ %s", cmd)
    proc = subprocess.run(
        cmd, shell=True, cwd=str(cwd), env=env,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    log.info(proc.stdout)
    if proc.returncode != 0:
        log.error("ERROR: '%s' exited with code %d.", label, proc.returncode)
        sys.exit(proc.returncode)
    return proc.stdout


def _check_mesh_quality(output: str) -> None:
    """
    Parse checkMesh output; abort if maxNonOrtho > 70 ° or maxSkewness > 4.
    These thresholds catch the pathological cells that crash pimpleFoam later.
    """
    bad = False

    m = re.search(r"Max non-orthogonality\s*=\s*([\d.]+)", output)
    if m and float(m.group(1)) > 70:
        print(
            f"ERROR: maxNonOrtho = {float(m.group(1)):.1f} ° > 70 ° threshold.",
            file=sys.stderr,
        )
        bad = True

    m = re.search(r"Max skewness\s*=\s*([\d.]+)", output)
    if m and float(m.group(1)) > 4:
        print(
            f"ERROR: maxSkewness = {float(m.group(1)):.2f} > 4 threshold.",
            file=sys.stderr,
        )
        bad = True

    if bad:
        log.error("Mesh quality check FAILED. Inspect the mesh in ParaView before proceeding.")
        sys.exit(1)

    log.info("[vortex-cfd] Mesh quality check PASSED.")


def _run_snappy(case_dir: Path, env: dict | None, inlet_params: dict) -> None:
    """
    Run snappyHexMesh, retrying up to 2 times with a progressively finer
    background mesh (target halved each time) if SHM exits non-zero.

    The retry re-renders only blockMeshDict and re-runs blockMesh before
    each SHM attempt — all other templates are left untouched.
    """
    from jinja2 import Environment, FileSystemLoader
    from .case_builder import _background_cell_counts

    target = 0.002   # 2 mm starting target
    bbox = inlet_params.get("_bbox")
    ctx = inlet_params.get("_ctx")
    tmpl_dir = inlet_params.get("_tmpl_dir")

    for attempt in range(3):
        cmd = "snappyHexMesh -overwrite"
        label = f"snappyHexMesh (attempt {attempt + 1}/3)"
        log.info("\n[vortex-cfd] %s", label)
        log.info("  $ %s", cmd)
        proc = subprocess.run(
            cmd, shell=True, cwd=str(case_dir), env=env,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
        )
        log.info(proc.stdout)
        if proc.returncode == 0:
            return

        if attempt == 2:
            log.error("snappyHexMesh failed after 3 attempts. Inspect the mesh setup.")
            sys.exit(proc.returncode)

        # Retry: halve target cell size, re-render blockMeshDict, re-run blockMesh
        target /= 2.0
        log.warning("snappyHexMesh failed (code %d). Retrying with target %.4f m ...",
                    proc.returncode, target)

        if bbox and ctx and tmpl_dir:
            counts = _background_cell_counts(bbox, target)
            new_ctx = {**ctx, "nx": counts["nx"], "ny": counts["ny"], "nz": counts["nz"]}
            jinja_env = Environment(
                loader=FileSystemLoader(tmpl_dir), keep_trailing_newline=True
            )
            rendered = jinja_env.get_template("system/blockMeshDict.j2").render(**new_ctx)
            (case_dir / "system" / "blockMeshDict").write_text(rendered, encoding="utf-8")

        _run("blockMesh", case_dir, env, "blockMesh (retry)")


def _read_inlet_geometry(case_dir: Path) -> tuple:
    """
    Read inlet patch face centres from the just-built mesh via pyvista.
    Returns (face_centers, ...) — same pattern as postprocess.read_wss_series.
    """
    import pyvista as pv
    foam_files = sorted(case_dir.glob("*.foam"))
    if not foam_files:
        raise FileNotFoundError(f"No .foam file in {case_dir}")
    reader = pv.OpenFOAMReader(str(foam_files[0]))
    reader.enable_all_patch_arrays()
    reader.cell_to_point_creation = False
    reader.set_active_time_value(reader.time_values[0])
    mb = reader.read()
    boundary = mb.get("boundary") or mb
    inlet = None
    for key in (boundary.keys() if hasattr(boundary, "keys") else []):
        if key == "inlet":
            inlet = boundary[key]
            break
    if inlet is None:
        raise KeyError("'inlet' patch not found in the mesh after snappyHexMesh")
    return np.array(inlet.cell_centers().points)


def run_pipeline(
    case_dir: Path,
    of_env: dict,
    cores: int | None,
    cycles: int | None = None,
    postprocess_metrics: bool = False,
    womersley: bool = False,
    inlet_params: dict | None = None,
) -> None:
    """
    Full Phase A/B pipeline:
      surfaceFeatureExtract
      blockMesh
      snappyHexMesh
      [write Womersley boundaryData if --womersley]
      checkMesh
      [decomposePar -force]
      pimpleFoam [-parallel]
      [reconstructPar]

    When postprocess_metrics is True the wallShearStress field is written live
    during the solve (function objects in controlDict) and Phase C metrics are
    computed after reconstructPar.

    When womersley is True the Womersley inlet profile data is written to
    constant/boundaryData/inlet/ after snappyHexMesh and before checkMesh.
    """
    if cores is None:
        cores = multiprocessing.cpu_count()

    env = _build_env(of_env)
    parallel = cores > 1
    mpi = f"mpirun -np {cores} " if parallel else ""

    _run("surfaceFeatureExtract",           case_dir, env, "surfaceFeatureExtract")
    _run("blockMesh",                        case_dir, env, "blockMesh")

    # snappyHexMesh runs in serial regardless of --cores. The parallel load-
    # balancing pass in v2406 has a segfault bug (fvMeshDistribute::repatch)
    # triggered when unbalance > 0.1 during shell refinement. Serial meshing
    # is safe and fast enough for typical ICA geometries (~2 min on 6 cores
    # would only save ~30 s). Parallel solving below is unaffected.
    _run_snappy(case_dir, env, inlet_params or {})

    # Womersley: read actual face positions from the meshed inlet patch and
    # write constant/boundaryData/inlet/<t>/U for timeVaryingMappedFixedValue.
    # Must happen after SHM (face positions only exist then) and before checkMesh.
    if womersley and inlet_params:
        from . import womersley as wom
        from .waveform import T_CYCLE
        log.info("\n[vortex-cfd] Writing Womersley boundary data ...")
        face_centers = _read_inlet_geometry(case_dir)
        t_abs, U_vecs = wom.womersley_velocities(
            face_centers=face_centers,
            centroid=inlet_params["centroid"],
            normal=inlet_params["normal"],
            radius=inlet_params["radius"],
            mean_velocity=inlet_params["mean_velocity"],
            waveform=inlet_params["waveform"],
            nu=inlet_params["nu"],
            t_cycle=T_CYCLE,
        )
        wom.write_boundary_data(case_dir, "inlet", t_abs, face_centers, U_vecs)
        log.info("  Written %d time steps to constant/boundaryData/inlet/", len(t_abs))

    check_out = _run("checkMesh",            case_dir, env, "checkMesh")
    _check_mesh_quality(check_out)

    if parallel:
        _run("decomposePar -force",          case_dir, env, "decomposePar (pimpleFoam)")

    pimple_cmd = f"{mpi}pimpleFoam" + (" -parallel" if parallel else "")
    _run(pimple_cmd,                         case_dir, env, "pimpleFoam")

    if parallel:
        _run("reconstructPar",               case_dir, env, "reconstructPar")

    if postprocess_metrics:
        postprocess.compute_metrics(case_dir, cycles=cycles)

    log.info("\n[vortex-cfd] Pipeline complete.")
    log.info("  Case : %s", case_dir)
    log.info("  Open : %s in ParaView", case_dir / (case_dir.name + ".foam"))


def generate_wss_fields(case_dir: Path, of_env: dict) -> None:
    """
    Generate the wallShearStress field from existing U snapshots without
    re-simulating, using pimpleFoam's -postProcess mode.  Used by the
    standalone path on a case that was solved before --postprocess existed.
    Operates on reconstructed time directories in the case root.
    """
    env = _build_env(of_env)
    _run(
        "pimpleFoam -postProcess -func wallShearStress",
        case_dir, env, "wallShearStress (-postProcess)",
    )


def run_postprocess_only(case_dir: Path, of_env: dict, cycles: int | None = None) -> None:
    """
    Standalone Phase C: compute biomarkers on an already-solved case.
    Generates the wallShearStress field first if it is not already present,
    then computes and writes metrics_report.json.
    """
    case_dir = Path(case_dir)

    # Has wallShearStress already been written into any snapshot?
    has_wss = any(case_dir.glob("[0-9]*/wallShearStress"))
    if not has_wss:
        log.info("[vortex-cfd] wallShearStress not found — generating via -postProcess.")
        generate_wss_fields(case_dir, of_env)

    postprocess.compute_metrics(case_dir, cycles=cycles)

    log.info("\n[vortex-cfd] Post-processing complete.")
    log.info("  Case   : %s", case_dir)
    log.info("  Report : %s", case_dir / "metrics_report.json")
