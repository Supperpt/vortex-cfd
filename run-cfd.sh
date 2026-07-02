#!/usr/bin/env bash
# Launcher: sources the newest accepted OpenFOAM ESI installation found, then
# runs the vortex-cfd CLI inside the vortex-aneurysm conda environment.
set -e

_sourced_foam=0
for ver in 2512 2506 2412 2406; do
    for prefix in /opt /usr/lib/openfoam; do
        bashrc="${prefix}/openfoam${ver}/etc/bashrc"
        if [ -f "$bashrc" ]; then
            # shellcheck disable=SC1090
            set +e; source "$bashrc" 2>/dev/null; set -e
            _sourced_foam=1
            break 2
        fi
    done
done

if [ "$_sourced_foam" -eq 0 ]; then
    echo "INFO: No OpenFOAM installation found in standard paths." >&2
    echo "      Continuing — env_check.py will abort with a clear message if needed." >&2
fi

# Locate python: prefer active conda env, then .venv, then conda run
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ -f "$SCRIPT_DIR/.venv/bin/python" ]; then
    exec "$SCRIPT_DIR/.venv/bin/python" -m vortex_cfd "$@"
elif [ "${CONDA_DEFAULT_ENV:-}" = "vortex-aneurysm" ] && [ -n "$CONDA_PREFIX" ]; then
    exec "$CONDA_PREFIX/bin/python" -m vortex_cfd "$@"
else
    # Active env is not vortex-aneurysm — locate the conda base and activate it.
    # Don't hardcode a single distribution: miniforge/mambaforge/anaconda all
    # live under different directories. Prefer `conda info --base` if conda is on
    # PATH, then fall back to the common install locations.
    CONDA_SH=""
    if command -v conda >/dev/null 2>&1; then
        _base="$(conda info --base 2>/dev/null)"
        [ -n "$_base" ] && [ -f "$_base/etc/profile.d/conda.sh" ] && CONDA_SH="$_base/etc/profile.d/conda.sh"
    fi
    if [ -z "$CONDA_SH" ]; then
        for _base in "$HOME/miniforge3" "$HOME/mambaforge" "$HOME/miniconda3" "$HOME/anaconda3" "/opt/conda"; do
            if [ -f "$_base/etc/profile.d/conda.sh" ]; then
                CONDA_SH="$_base/etc/profile.d/conda.sh"
                break
            fi
        done
    fi

    if [ -n "$CONDA_SH" ]; then
        source "$CONDA_SH"
        conda activate vortex-aneurysm
        exec "$CONDA_PREFIX/bin/python" -m vortex_cfd "$@"
    else
        echo "ERROR: Could not find conda or .venv. Activate vortex-aneurysm or run setup.sh first." >&2
        exit 1
    fi
fi
