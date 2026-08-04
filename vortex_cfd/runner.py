"""Orchestrate the full OpenFOAM meshing and solving pipeline."""
from __future__ import annotations

import multiprocessing
import os
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

from . import postprocess


def _inlet_face_centers(case_dir: Path, patch: str = "inlet"):
    """
    Face centres of the meshed inlet patch, read after snappyHexMesh.

    Uses the direct ``boundary`` lookup rather than
    ``postprocess._wall_block``: that helper's fallback path only recognises
    blocks carrying ``wallShearStress``, which does not exist before the solve.
    """
    import numpy as np
    import pyvista as pv

    foam_files = sorted(Path(case_dir).glob("*.foam"))
    if not foam_files:
        raise FileNotFoundError(f"No .foam file in {case_dir}")

    reader = pv.OpenFOAMReader(str(foam_files[0]))
    reader.enable_all_patch_arrays()
    reader.cell_to_point_creation = False
    mb = reader.read()

    boundary = mb["boundary"] if "boundary" in mb.keys() else mb
    if not hasattr(boundary, "keys") or patch not in boundary.keys():
        available = list(boundary.keys()) if hasattr(boundary, "keys") else []
        raise KeyError(
            f"patch '{patch}' not found in the mesh after snappyHexMesh "
            f"(available: {available})"
        )
    return np.asarray(boundary[patch].cell_centers().points, dtype=float)


def _write_womersley_boundary_data(case_dir: Path, inlet_params: dict | None) -> None:
    """
    Build the Womersley inlet profile on the meshed inlet faces and write
    constant/boundaryData/inlet/.

    Aborts the run on failure: unlike a post-processing diagnostic, a missing
    boundaryData directory would leave the solver with an inlet BC it cannot
    read, so failing here is far cheaper than failing hours into the solve.
    """
    from . import womersley as wom

    if not inlet_params:
        print("ERROR: --womersley requires inlet geometry from build_case, "
              "but none was supplied.", file=sys.stderr)
        _log_append("\n## Womersley boundaryData — FAILED (no inlet parameters)\n")
        sys.exit(1)

    print("\n[vortex-cfd] Writing Womersley inlet boundaryData")
    try:
        centers = _inlet_face_centers(case_dir)
        times, vectors = wom.womersley_velocities(
            face_centers=centers,
            centroid=inlet_params["centroid"],
            normal=inlet_params["normal"],
            radius=inlet_params["radius"],
            mean_velocity=inlet_params["mean_velocity"],
            waveform=inlet_params["waveform"],
            nu=inlet_params["nu"],
            cycles=inlet_params["cycles"],
        )
        wom.write_boundary_data(case_dir, "inlet", times, centers, vectors)
    except (FileNotFoundError, KeyError, ValueError) as e:
        print(f"ERROR: could not write Womersley boundaryData: {e}", file=sys.stderr)
        _log_append(f"\n## Womersley boundaryData — FAILED\n\n{e}\n")
        sys.exit(1)

    print(f"  {len(times)} time steps x {len(centers)} inlet faces "
          f"-> constant/boundaryData/inlet/")
    # A pure-Python step is not logged by _run(), so record it explicitly.
    _log_append(
        f"\n## Womersley boundaryData — OK\n\n"
        f"- Inlet faces: {len(centers)}\n"
        f"- Time steps: {len(times)} (0 to {times[-1]:.5f} s)\n"
        f"- Radius: {inlet_params['radius']:.6g} m, "
        f"mean velocity {inlet_params['mean_velocity']:.4g} m/s\n"
    )


# Path to the current run's markdown log, or None when logging is off.
# Set by _log_init() at the start of a pipeline; appended to by _run() and
# _check_mesh_quality() as the run proceeds (written incrementally so the log
# survives a mid-run abort).
_RUN_LOG: "Path | None" = None


def _log_init(case_dir: Path, pipeline: str) -> None:
    """Start a fresh run_log.md in the case directory."""
    global _RUN_LOG
    _RUN_LOG = Path(case_dir) / "run_log.md"
    try:
        _RUN_LOG.write_text(
            f"# vortex-cfd run log\n\n"
            f"- **Case:** `{case_dir}`\n"
            f"- **Started:** {datetime.now().isoformat(timespec='seconds')}\n"
            f"- **Pipeline:** {pipeline}\n"
        )
    except OSError:
        _RUN_LOG = None  # never let logging break a run


def _log_append(markdown: str) -> None:
    """Append a markdown fragment to the run log, if logging is active."""
    if _RUN_LOG is None:
        return
    try:
        with _RUN_LOG.open("a") as f:
            f.write(markdown)
    except OSError:
        pass


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
    print(f"\n[vortex-cfd] {label}")
    print(f"  $ {cmd}")
    proc = subprocess.run(
        cmd, shell=True, cwd=str(cwd), env=env,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    print(proc.stdout, end="")
    status = "OK" if proc.returncode == 0 else f"FAILED (exit {proc.returncode})"
    _log_append(
        f"\n## {label} — {status}\n\n"
        f"```sh\n{cmd}\n```\n\n"
        f"```\n{proc.stdout.rstrip()}\n```\n"
    )
    if proc.returncode != 0:
        print(
            f"\nERROR: '{label}' exited with code {proc.returncode}.",
            file=sys.stderr,
        )
        sys.exit(proc.returncode)
    return proc.stdout


def _check_mesh_quality(output: str) -> None:
    """
    Parse checkMesh output; abort if maxNonOrtho > 70 ° or maxSkewness > 20.
    These thresholds catch the pathological cells that crash pimpleFoam later.

    The skewness limit mirrors OpenFOAM's own standard: the meshing target is
    maxInternalSkewness 4, but boundary faces are allowed up to
    maxBoundarySkewness 20. checkMesh reports a single global "Max skewness" and
    flags anything > 4, so a handful of boundary-layer faces in the 4–20 range
    routinely trip an abort at 4 even though the mesh is usable (e.g. AA_011:
    max 6.66 on 29/3.12 M faces). 20 matches OpenFOAM's usable boundary limit;
    non-orthogonality and negative cell volumes remain the genuinely fatal checks.
    """
    bad = False

    # checkMesh prints "Mesh non-orthogonality Max: 69.96 average: 12.28"
    # (note the word order and ':' — not "Max non-orthogonality = ").
    ortho_m = re.search(r"non-orthogonality Max:\s*([\d.]+)", output)
    if ortho_m and float(ortho_m.group(1)) > 70:
        print(
            f"ERROR: maxNonOrtho = {float(ortho_m.group(1)):.1f} ° > 70 ° threshold.",
            file=sys.stderr,
        )
        bad = True

    skew_m = re.search(r"Max skewness\s*=\s*([\d.]+)", output)
    if skew_m and float(skew_m.group(1)) > 20:
        print(
            f"ERROR: maxSkewness = {float(skew_m.group(1)):.2f} > 20 threshold.",
            file=sys.stderr,
        )
        bad = True

    # Fail-safe: if neither metric was found, the output format changed (new
    # OpenFOAM release, locale, or a truncated log). Don't let the gate pass
    # silently — warn loudly so a bad mesh isn't waved through.
    if ortho_m is None and skew_m is None:
        print(
            "WARNING: could not parse non-orthogonality or skewness from checkMesh "
            "output — the mesh quality gate was NOT enforced. Inspect the mesh "
            "manually before trusting the solution.",
            file=sys.stderr,
        )
        _log_append("\n## Mesh quality gate — NOT ENFORCED\n\n"
                    "checkMesh output could not be parsed (format change?).\n")

    if bad:
        print(
            "\nMesh quality check FAILED. Inspect the mesh in ParaView before proceeding.",
            file=sys.stderr,
        )
        _log_append("\n## Mesh quality gate — FAILED\n\n"
                    "Aborted before solving (see checkMesh output above).\n")
        sys.exit(1)

    print("[vortex-cfd] Mesh quality check PASSED.")
    _log_append("\n## Mesh quality gate — PASSED\n")


def run_pipeline(
    case_dir: Path,
    of_env: dict,
    cores: int | None,
    cycles: int | None = None,
    postprocess_metrics: bool = False,
    neck_metrics: bool = False,
    womersley: bool = False,
    inlet_params: dict | None = None,
) -> None:
    """
    Full Phase A pipeline:
      surfaceFeatureExtract
      blockMesh
      [decomposePar]
      snappyHexMesh [-parallel]
      [reconstructParMesh -constant]
      checkMesh
      [decomposePar -force]
      pimpleFoam [-parallel]
      [reconstructPar]

    When postprocess_metrics is True the wallShearStress field is written live
    during the solve (function objects in controlDict) and Phase C metrics are
    computed after reconstructPar.
    """
    if cores is None:
        cores = multiprocessing.cpu_count()

    env = _build_env(of_env)
    parallel = cores > 1
    mpi = f"mpirun -np {cores} " if parallel else ""

    _log_init(case_dir, f"full pipeline (cores={cores}, cycles={cycles}, "
                        f"postprocess={postprocess_metrics}, womersley={womersley})")

    _run("surfaceFeatureExtract",           case_dir, env, "surfaceFeatureExtract")
    _run("blockMesh",                        case_dir, env, "blockMesh")

    # snappyHexMesh runs in serial regardless of --cores. The parallel load-
    # balancing pass in v2406 has a segfault bug (fvMeshDistribute::repatch)
    # triggered when unbalance > 0.1 during shell refinement. Serial meshing
    # is safe and fast enough for typical ICA geometries (~2 min on 6 cores
    # would only save ~30 s). Parallel solving below is unaffected.
    _run("snappyHexMesh -overwrite",        case_dir, env, "snappyHexMesh")

    # Womersley boundaryData must be written here: the inlet face centres only
    # exist once the mesh does, and the files must be in place before
    # decomposePar distributes the case below.
    if womersley:
        _write_womersley_boundary_data(case_dir, inlet_params)

    check_out = _run("checkMesh",            case_dir, env, "checkMesh")
    _check_mesh_quality(check_out)

    if parallel:
        _run("decomposePar -force",          case_dir, env, "decomposePar (pimpleFoam)")

    pimple_cmd = f"{mpi}pimpleFoam" + (" -parallel" if parallel else "")
    _run(pimple_cmd,                         case_dir, env, "pimpleFoam")

    if parallel:
        _run("reconstructPar",               case_dir, env, "reconstructPar")

    if postprocess_metrics:
        postprocess.compute_metrics(case_dir, cycles=cycles, neck_metrics=neck_metrics)
        _log_append("\n## Post-processing — metrics_report.json written\n")
        if neck_metrics:
            _log_append("Neck inflow metrics computed (experimental, unvalidated).\n")

    _log_append(f"\n## Pipeline complete — {datetime.now().isoformat(timespec='seconds')}\n")
    print(f"\n[vortex-cfd] Pipeline complete.")
    print(f"  Case : {case_dir}")
    print(f"  Open : {case_dir / (case_dir.name + '.foam')} in ParaView")
    print(f"  Log  : {case_dir / 'run_log.md'}")


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


def run_postprocess_only(case_dir: Path, of_env: dict, cycles: int | None = None,
                         neck_metrics: bool = False) -> None:
    """
    Standalone Phase C: compute biomarkers on an already-solved case.
    Generates the wallShearStress field first if it is not already present,
    then computes and writes metrics_report.json.
    """
    case_dir = Path(case_dir)
    _log_init(case_dir, "post-process only")

    # Has wallShearStress already been written into any snapshot?
    has_wss = any(case_dir.glob("[0-9]*/wallShearStress"))
    if not has_wss:
        print("[vortex-cfd] wallShearStress not found — generating via -postProcess.")
        generate_wss_fields(case_dir, of_env)

    postprocess.compute_metrics(case_dir, cycles=cycles, neck_metrics=neck_metrics)

    print(f"\n[vortex-cfd] Post-processing complete.")
    print(f"  Case   : {case_dir}")
    print(f"  Report : {case_dir / 'metrics_report.json'}")
