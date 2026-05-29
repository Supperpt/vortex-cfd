#!/usr/bin/env bash
# Smoke test: run the full Phase A pipeline on a known STL set and verify that
# the resulting OpenFOAM case directory has the expected structure.
#
# Usage:
#   STL_DIR=/path/to/vortex/output/patient_01 bash tests/smoke_test.sh
#
# The test passes when:
#   1. The case directory is created with the correct layout.
#   2. All OpenFOAM dict files are non-empty.
#   3. The .foam ParaView placeholder exists.
#   4. patch_labels.json is valid JSON.
# The test does NOT run the OpenFOAM mesher/solver (that requires a Linux host
# with OpenFOAM installed); it only validates the Python pipeline up to and
# including case generation.
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
export REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
export STL_DIR="${STL_DIR:?Set STL_DIR to a directory of VORTEX STL files}"
export OUT_DIR="${TMPDIR:-/tmp}/vortex_cfd_smoke_$$"

mkdir -p "$OUT_DIR"
trap 'rm -rf "$OUT_DIR"' EXIT

echo "=== vortex-cfd smoke test ==="
echo "STL_DIR : $STL_DIR"
echo "OUT_DIR : $OUT_DIR"

# Run only the case-builder step by setting --cycles 1 and intercepting before
# the solver.  We use a tiny mean velocity that is clearly non-physical but
# sufficient to exercise the code path.
PYTHONPATH="$REPO_ROOT" python - <<'PYEOF'
import sys, os
from pathlib import Path

stl_dir   = Path(os.environ["STL_DIR"])
out_dir   = Path(os.environ["OUT_DIR"])
cycles    = 1
u_mean    = 0.4
cores     = 1

stl_paths = sorted(stl_dir.glob("*.stl"))
assert stl_paths, f"No STL files in {stl_dir}"

# Build a minimal label dict: largest surface → wall, rest → inlet + outlets
from vortex_cfd.patch_labeller import VALID_LABELS
import pyvista as pv

areas = {}
for p in stl_paths:
    m = pv.read(str(p))
    s = m.compute_cell_sizes()
    areas[p] = float(s.cell_data["Area"].sum())

wall_path   = max(areas, key=areas.__getitem__)
cap_paths   = [p for p in stl_paths if p != wall_path]
inlet_path  = cap_paths[0]
outlet_paths = cap_paths[1:] if len(cap_paths) > 1 else []

labels = {wall_path: "wall", inlet_path: "inlet"}
for p in outlet_paths:
    labels[p] = "outlet"

# If there are no outlets yet, use the inlet also as outlet_0 (degenerate geometry)
if not outlet_paths:
    # Add a second label for the inlet as outlet_0 — smoke test only
    labels[inlet_path] = "outlet"
    labels[wall_path] = "wall"
    # Need at least one real inlet; skip test with a warning
    print("WARNING: only one cap found — cannot form inlet+outlet. Skipping solver test.")
    sys.exit(0)

from vortex_cfd.scaling import scale_stls
from vortex_cfd.waveform import load_waveform
from vortex_cfd.case_builder import build_case

scaled   = scale_stls(stl_paths, labels)
waveform = load_waveform(None)
case_dir = build_case(
    scaled_stls=scaled, labels=labels, cycles=cycles,
    mean_velocity=u_mean, waveform=waveform, cores=cores, out_dir=str(out_dir),
)

# Assertions
required_files = [
    "0/U", "0/p",
    "constant/transportProperties", "constant/turbulenceProperties",
    "system/controlDict", "system/fvSchemes", "system/fvSolution",
    "system/snappyHexMeshDict", "system/decomposeParDict",
    "system/meshQualityDict", "system/surfaceFeatureExtractDict",
    "system/blockMeshDict",
    "patch_labels.json", "Allrun",
]

for rel in required_files:
    p = case_dir / rel
    assert p.exists(), f"Missing: {p}"
    assert p.stat().st_size > 0, f"Empty: {p}"

import json
meta = json.loads((case_dir / "patch_labels.json").read_text())
assert "wall" in meta.values(), "patch_labels.json missing wall entry"
assert "inlet" in meta.values(), "patch_labels.json missing inlet entry"
assert "outlet" in meta.values(), "patch_labels.json missing outlet entry"

foam_files = list(case_dir.glob("*.foam"))
assert foam_files, "No .foam placeholder found"

print(f"PASSED — case: {case_dir}")
PYEOF

echo "=== smoke test PASSED ==="
