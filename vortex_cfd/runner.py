"""Orchestrate the full OpenFOAM meshing and solving pipeline."""

import multiprocessing
import os
import re
import subprocess
import sys
from pathlib import Path


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
    if proc.returncode != 0:
        print(
            f"\nERROR: '{label}' exited with code {proc.returncode}.",
            file=sys.stderr,
        )
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
        print(
            "\nMesh quality check FAILED. Inspect the mesh in ParaView before proceeding.",
            file=sys.stderr,
        )
        sys.exit(1)

    print("[vortex-cfd] Mesh quality check PASSED.")


def run_pipeline(case_dir: Path, of_env: dict, cores: int | None) -> None:
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
    """
    if cores is None:
        cores = multiprocessing.cpu_count()

    env = _build_env(of_env)
    parallel = cores > 1
    mpi = f"mpirun -np {cores} " if parallel else ""

    _run("surfaceFeatureExtract",           case_dir, env, "surfaceFeatureExtract")
    _run("blockMesh",                        case_dir, env, "blockMesh")

    if parallel:
        _run("decomposePar",                 case_dir, env, "decomposePar (snappyHexMesh)")

    snap_cmd = f"{mpi}snappyHexMesh -overwrite" + (" -parallel" if parallel else "")
    _run(snap_cmd,                           case_dir, env, "snappyHexMesh")

    if parallel:
        _run("reconstructParMesh -constant", case_dir, env, "reconstructParMesh")

    check_out = _run("checkMesh",            case_dir, env, "checkMesh")
    _check_mesh_quality(check_out)

    if parallel:
        _run("decomposePar -force",          case_dir, env, "decomposePar (pimpleFoam)")

    pimple_cmd = f"{mpi}pimpleFoam" + (" -parallel" if parallel else "")
    _run(pimple_cmd,                         case_dir, env, "pimpleFoam")

    if parallel:
        _run("reconstructPar",               case_dir, env, "reconstructPar")

    print(f"\n[vortex-cfd] Pipeline complete.")
    print(f"  Case : {case_dir}")
    print(f"  Open : {case_dir / (case_dir.name + '.foam')} in ParaView")
