"""Detect and validate an active OpenFOAM ESI (openfoam.com) environment."""

import os
import shutil
import subprocess
import sys
from pathlib import Path

ACCEPTED_VERSIONS = {"2406", "2412", "2506", "2512"}
_SEARCH_ROOTS = [
    "/opt",
    "/usr/lib/openfoam",
]


def _active_version() -> str | None:
    """Return WM_PROJECT_VERSION if blockMesh is already on PATH."""
    if shutil.which("blockMesh") is None:
        return None
    return os.environ.get("WM_PROJECT_VERSION")


def _source_and_capture(bashrc: Path) -> dict[str, str]:
    result = subprocess.run(
        f'bash -c "source {bashrc} && env"',
        shell=True, capture_output=True, text=True,
    )
    if result.returncode != 0:
        return {}
    env: dict[str, str] = {}
    for line in result.stdout.splitlines():
        k, _, v = line.partition("=")
        if k:
            env[k] = v
    return env


def _normalise(version: str) -> str:
    """Strip leading 'v' so '2512' and 'v2512' both compare against ACCEPTED_VERSIONS."""
    return version.lstrip("v")


def check_openfoam() -> dict:
    """
    Return {'version', 'root', 'env', 'bashrc'} for the active OpenFOAM ESI install.
    'env' is empty when OpenFOAM is already sourced in the calling shell.
    Exits with a helpful message when no acceptable installation is found.
    """
    version = _active_version()
    if version:
        norm = _normalise(version)
        if norm not in ACCEPTED_VERSIONS:
            print(
                f"WARNING: OpenFOAM version '{version}' is outside the tested set "
                f"({', '.join(sorted(ACCEPTED_VERSIONS))}). Proceeding anyway.",
                file=sys.stderr,
            )
        root = os.environ.get("WM_PROJECT_DIR", os.environ.get("FOAM_INST_DIR", ""))
        return {"version": norm, "root": root, "env": {}, "bashrc": None}

    # OpenFOAM is not sourced — search standard installation paths.
    for root_dir in _SEARCH_ROOTS:
        for ver in sorted(ACCEPTED_VERSIONS, reverse=True):
            bashrc = Path(root_dir) / f"openfoam{ver}" / "etc" / "bashrc"
            if bashrc.exists():
                env = _source_and_capture(bashrc)
                found_ver = _normalise(env.get("WM_PROJECT_VERSION", ""))
                if found_ver in ACCEPTED_VERSIONS:
                    print(
                        f"NOTE: OpenFOAM {found_ver} found at {bashrc.parent.parent} "
                        f"but not sourced in the current shell.\n"
                        f"      For manual inspection: source {bashrc}",
                        file=sys.stderr,
                    )
                    return {
                        "version": found_ver,
                        "root": str(bashrc.parent.parent),
                        "env": env,
                        "bashrc": str(bashrc),
                    }

    print(
        "ERROR: No OpenFOAM ESI installation found.\n"
        "  Accepted versions : " + ", ".join(sorted(ACCEPTED_VERSIONS)) + "\n"
        "  Searched under    : " + ", ".join(_SEARCH_ROOTS) + "\n"
        "  Install from https://openfoam.com and source its etc/bashrc before running.",
        file=sys.stderr,
    )
    sys.exit(1)
