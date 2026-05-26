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
            source "$bashrc"
            _sourced_foam=1
            break 2
        fi
    done
done

if [ "$_sourced_foam" -eq 0 ]; then
    echo "INFO: No OpenFOAM installation found in standard paths." >&2
    echo "      Continuing — env_check.py will abort with a clear message if needed." >&2
fi

exec conda run -n vortex-aneurysm python -m vortex_cfd.cli "$@"
