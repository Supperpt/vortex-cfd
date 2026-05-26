#!/usr/bin/env bash
# Add vortex-cfd Python deps to the shared vortex-aneurysm conda environment.
set -e
conda run -n vortex-aneurysm pip install -r requirements.txt
echo "Done. Activate with: conda activate vortex-aneurysm"
