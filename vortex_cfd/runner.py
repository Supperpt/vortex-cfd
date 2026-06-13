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
    m = re.search(r"non-orthogonality Max:\s*([\d.]+)", output)
    if m and float(m.group(1)) > 70:
        print(
            f"ERROR: maxNonOrtho = {float(m.group(1)):.1f} ° > 70 ° threshold.",
            file=sys.stderr,
        )
        bad = True

    m = re.search(r"Max skewness\s*=\s*([\d.]+)", output)
    if m and float(m.group(1)) > 20:
        print(
            f"ERROR: maxSkewness = {float(m.group(1)):.2f} > 20 threshold.",
            file=sys.stderr,
        )
        bad = True

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
                        f"postprocess={postprocess_metrics})")

    _run("surfaceFeatureExtract",           case_dir, env, "surfaceFeatureExtract")
    _run("blockMesh",                        case_dir, env, "blockMesh")

    # snappyHexMesh runs in serial regardless of --cores. The parallel load-
    # balancing pass in v2406 has a segfault bug (fvMeshDistribute::repatch)
    # triggered when unbalance > 0.1 during shell refinement. Serial meshing
    # is safe and fast enough for typical ICA geometries (~2 min on 6 cores
    # would only save ~30 s). Parallel solving below is unaffected.
    _run("snappyHexMesh -overwrite",        case_dir, env, "snappyHexMesh")

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
        _log_append("\n## Post-processing — metrics_report.json written\n")

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


def run_postprocess_only(case_dir: Path, of_env: dict, cycles: int | None = None) -> None:
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

    postprocess.compute_metrics(case_dir, cycles=cycles)

    print(f"\n[vortex-cfd] Post-processing complete.")
    print(f"  Case   : {case_dir}")
    print(f"  Report : {case_dir / 'metrics_report.json'}")
