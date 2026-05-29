# LLM.md — vortex-cfd session context

This file is the persistent handoff document between LLM sessions.
It is **not** user documentation — it records what has been built, why decisions were made, what is planned next, and every bug found with its resolution.

Update this file at the end of every session before closing.

---

## 1. Project overview

**One-sentence summary:** Automated end-to-end CFD pipeline that takes a patient STL from VORTEX and produces a complete, runnable OpenFOAM case with pulsatile Navier–Stokes solved and wall-shear-stress-ready mesh.

**PhD workflow position:**
```
DICOM (Angio-CT)  →[VORTEX]→  Watertight STL  →[vortex-cfd]→  WSS / OSI / TAWSS report
```

**Target platform:** Linux only. OpenFOAM ESI (openfoam.com), accepted versions: v2406, v2412, v2506, v2512. Python 3.10+ in the `vortex-aneurysm` conda environment (shared with VORTEX).

**VORTEX output contract:** A directory of STL files produced by `--split-patches`: one file for the lumen wall surface and one per capped opening. VMTK assigns cap IDs geometrically (no semantic labels); the user identifies inlet vs. outlets interactively at runtime.

---

## 2. Numerical choices (locked — do not change without updating README)

| Property | Value | Reason |
|---|---|---|
| Density | 1060 kg/m³ | Whole blood |
| Viscosity model | Newtonian, ν = 3.3 × 10⁻⁶ m²/s | Baseline; Carreau deferred to Phase D |
| Turbulence | Laminar | Re ≈ 200–400 in ICA parent artery |
| Solver | pimpleFoam | Transient, incompressible, PIMPLE |
| Time scheme | backward (2nd-order) | |
| Advection scheme | Gauss linearUpwind | Stable, ~2nd-order |
| Adaptive Δt | maxCo = 0.8 | Mandatory — peak systole is ~10× diastole |
| Cardiac period | T = 0.857 s | 70 bpm default |
| Inlet BC | flowRateInletVelocity | Uniform parabolic; Womersley in Phase B |
| Outlet BC | inletOutlet (U), fixedValue 0 (p) | Prevents recirculation instability |
| Wall BC | noSlip (U), zeroGradient (p) | Standard rigid-wall |
| BL layers | 4 prismatic, expansion 1.3, finalLayerThickness 0.3 | Required for WSS accuracy |
| BL refinement | level (3 4) → ~0.125 mm at wall | From 2 mm background cells |
| checkMesh abort | maxNonOrtho > 70°, maxSkewness > 4 | Cells that crash pimpleFoam later |

---

## 3. Test suite

Run with: `pytest` (from the repo root, after `pip install -e .`).

| File | What it tests | Tests |
|---|---|---|
| `tests/test_waveform.py` | Default waveform shape/normalisation, CSV loading | 15 |
| `tests/test_scaling.py` | mm detection, scaling factor, canonical name assignment | 14 |
| `tests/test_case_builder.py` | bbox+buffer, cell counts, inlet area, locationInMesh, waveform table, full case generation | 59 |
| `tests/test_env_check.py` | `_normalise`, accepted versions, `_active_version` with monkeypatching | 11 |

Total: **99 tests, 0 failures** (validated 2026-05-26 on Windows with synthetic STL geometry).

The tests require only Python + numpy + pyvista + jinja2 + pytest — no OpenFOAM installation needed.

---

## 4. Repository layout (current state)

```
vortex-cfd/
├── LLM.md                         ← this file
├── README.md                      ← user-facing documentation
├── pyproject.toml                 ← package definition + pytest config
├── requirements.txt               ← jinja2, numpy, pyvista, click, pytest
├── setup.sh                       ← pip install into vortex-aneurysm conda env
├── run-cfd.sh                     ← sources OpenFOAM, then runs the CLI
├── vortex_cfd/
│   ├── __init__.py
│   ├── __main__.py
│   ├── cli.py                     ← Click entry point
│   ├── env_check.py               ← detect/source OpenFOAM ESI
│   ├── patch_labeller.py          ← interactive wall/inlet/outlet labelling
│   ├── scaling.py                 ← mm→m detection and STL rewriting
│   ├── waveform.py                ← default ICA waveform + user CSV loader
│   ├── case_builder.py            ← geometry analysis + Jinja2 rendering
│   ├── runner.py                  ← pipeline orchestration (subprocess calls)
│   └── templates/

│       ├── 0/U.j2
│       ├── 0/p.j2
│       ├── Allrun.j2
│       ├── constant/
│       │   ├── transportProperties.j2
│       │   └── turbulenceProperties.j2
│       └── system/
│           ├── blockMeshDict.j2
│           ├── controlDict.j2
│           ├── decomposeParDict.j2
│           ├── fvSchemes.j2
│           ├── fvSolution.j2
│           ├── meshQualityDict.j2
│           ├── snappyHexMeshDict.j2
│           └── surfaceFeatureExtractDict.j2
└── tests/
    ├── conftest.py                ← pyvista-based synthetic STL fixtures
    ├── test_waveform.py           ← 15 tests
    ├── test_scaling.py            ← 14 tests
    ├── test_case_builder.py       ← 59 tests
    ├── test_env_check.py          ← 11 tests
    └── smoke_test.sh              ← shell-based end-to-end test (Linux only)
```

---

## 4. Implementation phases

### Phase A — MVP (COMPLETE as of 2026-05-26)

**Goal:** STLs in → runnable OpenFOAM case out → opens in ParaView.

**What was built:**
- Full Python package with Click CLI (`--stl-dir`, `--cycles`, `--mean-velocity`, `--waveform`, `--cores`, `--out-dir`)
- `env_check.py`: detects active OpenFOAM or searches standard paths; sources `etc/bashrc` and captures the env if needed
- `patch_labeller.py`: interactive prompt showing geometric description (centre, area, bounding box) for each STL; validates exactly 1 wall + 1 inlet + ≥1 outlet
- `scaling.py`: detects mm coordinates (bounding-box max > 1.0) and scales ×0.001 before saving to temp dir; writes canonical names (wall.stl, inlet.stl, outlet_0.stl, …)
- `waveform.py`: 3-harmonic Fourier analytical ICA waveform, normalised mean = 1; also loads user 2-column CSV
- `case_builder.py`: computes bounding box + 20% buffer, background cell counts (target 2 mm), `locationInMesh` (centre of wall STL bounding box), inlet area (pyvista), Q(t) table; renders all 13 templates; writes `patch_labels.json` AFTER the directory exists; touches `.foam` placeholder
- `runner.py`: orchestrates `surfaceFeatureExtract → blockMesh → [decomposePar] → snappyHexMesh [-parallel] → [reconstructParMesh] → checkMesh → [decomposePar -force] → pimpleFoam [-parallel] → [reconstructPar]`; aborts with clear message if checkMesh fails quality thresholds
- 13 OpenFOAM templates (Jinja2, validated to render without errors)
- `Allrun` shell script generated in each case directory for manual re-runs

**Not yet tested:** Against a real VORTEX STL on a Linux host with OpenFOAM installed. The `smoke_test.sh` tests only the Python case-generation step (no mesher/solver required).

**Success criterion (README):** A real patient STL goes in, a `case_XXX.foam` comes out that opens in ParaView and shows reasonable velocity fields.

---

### Phase B — Robustness and Womersley (PLANNED)

**Goal:** Production-grade inlet physics and hardened meshing.

Planned work items (in priority order):

1. **Womersley inlet profile** — `codedFixedValue` in `0/U` with Fourier decomposition of the waveform and Bessel functions for the radial profile. Required for academic publication. Needs a `--womersley` opt-in flag.

2. **Physiological validation** — Before meshing, check that `mean_velocity × inlet_area` gives a plausible cardiac output; warn (not abort) if outside typical ICA range (1–10 mL/s).

3. **Robust `locationInMesh`** — Current implementation uses the centre of the wall STL bounding box. This fails for highly curved or C-shaped vessels (centre of bbox is outside the lumen). Replace with: read inlet STL face centroids, compute the inlet face normal from area-weighted average, step one inlet-radius inward from the inlet centroid along that normal.

4. **snappyHexMesh retry logic** — If SHM exits non-zero (often due to insufficient background cells or an awkward geometry), automatically retry with a finer background mesh (halve the target cell size) up to 2 retries.

5. **Structured logging** — Replace `print()` calls with Python `logging` module; write a `vortex_cfd.log` inside the case directory.

6. **`--dry-run` flag** — Build the case directory but skip the OpenFOAM solver steps; useful for inspecting the mesh setup without waiting for a simulation.

---

### Phase C — Post-processing and reports (PLANNED)

**Goal:** WSS, TAWSS, OSI biomarkers as output.

Planned work items:

1. **OpenFOAM function objects** — Add `fieldAverage` and `wallShearStress` to `controlDict` so that `pimpleFoam` computes and writes `wallShearStressMean` and `UMean` during the last cycle.

2. **OSI computation** — OSI = 0.5 × (1 − |TAWSS| / ∫|WSS|dt) requires time-series WSS data. Options: (a) post-process with `foamCalc` after the run; (b) use a custom function object; (c) post-process with pyvista in Python. Option (c) is most portable.

3. **JSON report** — `metrics_report.json` with: mean TAWSS, max TAWSS, mean OSI, area fraction with OSI > 0.3, area fraction with TAWSS < 0.4 Pa (low-and-oscillatory risk zone), simulation metadata.

4. **Optional ParaView screenshots** — `pvbatch` script (generated by a template) to produce wall-coloured renders of WSS, OSI, TAWSS. Guarded by `--screenshots` flag.

5. **Validation ranges** (from prior iteration, `previous_cfd_iteration.md`):
   - WSS: 0–50 Pa
   - OSI: 0–0.5
   - Velocity: 0.1–1.0 m/s
   - Pressure: 0–200 Pa (relative to outlet)

---

### Phase D — Extensions (PLANNED)

1. **Carreau non-Newtonian viscosity** — Cho & Kensey (1991): μ₀ = 0.056, μ∞ = 0.0035 Pa·s, λ = 3.313 s, n = 0.3568. Opt-in via `--carreau` flag. Requires modifying `transportProperties.j2` and adding `CarreauYasuda` model.

2. **k-ω SST turbulence model** — Opt-in via `--turbulence` flag. Requires new `0/k.j2`, `0/omega.j2` templates and updated `turbulenceProperties.j2`.

3. **Mesh-independence helper** — Runs three mesh densities (coarse/medium/fine) and reports WSS convergence. Scheduled for after Phase C.

4. **FSI placeholder** — Reserved folder hooks for future fluid–structure interaction work.

---

## 5. Known issues and bugs

*No bugs recorded yet — Phase A just completed.*

**Format for future entries:**

```
### BUG-001 — [Short title]
- **Discovered:** YYYY-MM-DD, during [phase/activity]
- **Symptom:** [What went wrong / what error appeared]
- **Root cause:** [Why it happened]
- **Fix:** [What was changed, with file:line references]
- **Status:** FIXED / OPEN / DEFERRED
```

---

## 6. Design decisions log

Records non-obvious choices so future sessions don't re-litigate them.

### D-001 — mm→m scaling written back to STL before any other use
*2026-05-26*
The prior iteration (VesselForge_AutoCFD) applied scaling to the geometry in memory but forgot to rewrite patch-coordinate metadata files. This silently mapped the inlet patch to a point 1000× away from the mesh, producing runs that "completed" against the wrong geometry. In this implementation, `scaling.py` rewrites scaled copies to a temp directory immediately; everything downstream (case builder, SHM dict, patch labels) reads only the scaled copies.

### D-002 — `patch_labels.json` written after case directory exists
*2026-05-26*
Prior iteration wrote this file before the directory was created, leading to a race condition where the file was created in the current working directory instead. `case_builder.py` now writes it as the last step before returning.

### D-003 — `locationInMesh` = centre of wall STL bounding box
*2026-05-26*
Chosen for simplicity in Phase A. Works for typical ICA aneurysm geometries (roughly convex lumen). Known failure mode: highly curved vessels where the bounding-box centre falls inside the wall material. Phase B will replace this with an inlet-normal displacement method (see Phase B item 3).

### D-004 — `flowRateInletVelocity` with `outOfBounds repeat`
*2026-05-26*
Only one cardiac cycle is stored in the waveform table. OpenFOAM repeats it via `outOfBounds repeat`. This avoids duplicating the table N times for N cycles, and allows the cycle count to be changed without regenerating the case.

### D-005 — `pimpleFoam` with `nOuterCorrectors 2`
*2026-05-26*
The PIMPLE algorithm with 2 outer correctors gives a good balance between stability and cost at Co < 1. Increasing to 3 would be safer for very coarse meshes but adds 50% cost per timestep.

---

## 7. How to start a new session

1. Read this file (`LLM.md`) first.
2. Read `README.md` for user-facing requirements and numerical choices.
3. Check what changed since the last session: `git log --oneline -10`.
4. The current implementation phase is listed in Section 4.
5. Any open bugs are in Section 5.
6. Pick the first planned item from the next phase and implement it.

---

## ⚡ NEXT ACTION (start here)

**Run the full OpenFOAM pipeline against the real patient STL (Phase A validation).**

Everything is confirmed ready:
- OpenFOAM v2512 found at `/usr/lib/openfoam/openfoam2512/etc/bashrc` (picked up automatically by `run-cfd.sh`)
- STL files confirmed at `/home/diogovluz/Documentos/Programas/Programa VMTK/`
  - `output_wall.stl` → wall (31346 cells, area 598 mm²)
  - `input_cap_2.stl` → inlet (48 cells, area 6.5 mm²)
  - `output_cap_3.stl` → outlet (48 cells, area 1.8 mm²)
  - `output_cap_4.stl` → outlet (48 cells, area 3.2 mm²)
- Geometry is mm-scale (max extent ~36 mm) — auto-scaling ×0.001 will trigger correctly
- `vortex-aneurysm` conda env has all Python deps installed

**Command to run (single line, paste into terminal, interactive labelling):**

```bash
bash /home/diogovluz/Documentos/Programas/vortex-cfd/run-cfd.sh --stl-dir '/home/diogovluz/Documentos/Programas/Programa VMTK' --cycles 1 --mean-velocity 0.4 --cores 4 --out-dir /tmp/vortex_test
```

Use `--cycles 1` for the first test run (faster). When the interactive labeller prompts, assign:
- `output_wall.stl` → `wall`
- `input_cap_2.stl` → `inlet`
- `output_cap_3.stl` → `outlet`
- `output_cap_4.stl` → `outlet`

**Success criterion:** `case_XXXXXXXX_XXXXXX.foam` appears in `/tmp/vortex_test/` and opens in ParaView showing velocity fields. If anything fails, record it as a bug in Section 5 and fix before starting Phase B.

---

## 8. Session history

| Date | Session summary |
|---|---|
| 2026-05-26 | Phase A implementation: all Python modules + 13 OpenFOAM templates + Allrun script. Syntax-validated. Not yet run against real STLs on Linux. |
| 2026-05-26 | pytest suite: 99 tests across 4 files (waveform, scaling, case_builder, env_check). All pass on Windows with synthetic pyvista geometry. Added pyproject.toml. |
| 2026-05-29 | Phase A validation (partial): pytest confirmed 99/99 pass on Linux in vortex-aneurysm env. Fixed smoke_test.sh (OUT_DIR/STL_DIR/REPO_ROOT not exported — Python subprocess couldn't read them via os.environ). Smoke test PASSED with real VMTK STLs. OpenFOAM v2512 confirmed at standard path. Full mesher+solver run not yet executed — see NEXT ACTION above. |
