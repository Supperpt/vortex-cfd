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
elif [ -n "$CONDA_PREFIX" ]; then
    exec "$CONDA_PREFIX/bin/python" -m vortex_cfd "$@"
else
    # Try to find conda and use it
    CONDA_SH="$HOME/miniconda3/etc/profile.d/conda.sh"
    if [ -f "$CONDA_SH" ]; then
        source "$CONDA_SH"
        exec conda run -n vortex-aneurysm python -m vortex_cfd.cli "$@"
    else
        echo "ERROR: Could not find conda or .venv. Run setup first." >&2
        exit 1
    fi
fi
