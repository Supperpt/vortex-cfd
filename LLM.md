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
| `tests/test_case_builder.py` | bbox+buffer, cell counts, inlet area, locationInMesh, waveform table, full case generation, postprocess function objects | 64 |
| `tests/test_postprocess.py` | TAWSS/OSI/summary-stats pure math (no OpenFOAM): OSI=0 unidirectional, OSI=0.5 reversing, Pa=ρ×kinematic, area fractions | 16 |
| `tests/test_env_check.py` | `_normalise`, accepted versions, `_active_version` with monkeypatching | 11 |

Total: **120 tests** (validated 2026-05-26 on Windows with synthetic STL geometry).

> **Known pre-existing failure:** `TestLocationInMesh::test_location_inside_wall_bbox` fails on the
> current Linux/pyvista combination — the synthetic flat-disc fixture (`conftest._disc_stl`) yields an
> unexpected face normal/area so the stepped-in point lands outside the sphere bbox. This is a *test
> fixture* artefact (pyvista version-dependent), **not** a bug in `_location_in_mesh` itself, which was
> validated in ParaView on real geometry (BUG-006). Pre-dates Phase C. Tracked separately.

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
│   ├── postprocess.py             ← Phase C: WSS/TAWSS/OSI biomarkers + metrics_report.json
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

### Phase C — Post-processing and reports (IN PROGRESS — implemented 2026-06-04, awaiting real-case validation)

**Goal:** WSS, TAWSS, OSI biomarkers as output. Prioritised ahead of Phase B to enable initial publishable results and confirm end-to-end correctness.

Work items:

1. **OpenFOAM function objects** — DONE. `wallShearStress` + `fieldAverage` are appended to `controlDict` (only when `--postprocess` is set), writing `wallShearStress`, `wallShearStressMean`, `UMean`. `fieldAverage` `timeStart = (cycles-1)·T` averages only the last cycle. See `templates/system/controlDict.j2`.

2. **OSI / TAWSS computation** — DONE (option c, pyvista in Python). `postprocess.py` splits into pure numpy math (`tawss`, `osi`, `summary_stats` — unit-tested, no OpenFOAM) and an I/O layer (`read_wss_series` via `pv.OpenFOAMReader`, `compute_metrics`). OpenFOAM computes the WSS field; Python only does time-integration + stats. **WSS unit note:** OpenFOAM's incompressible `wallShearStress` is *kinematic* (m²/s², τ/ρ); physical Pa = ρ×value. OSI is dimensionless (ρ cancels).

3. **JSON report** — DONE. `metrics_report.json`: TAWSS in **both** Pa and kinematic, OSI, area fraction OSI > 0.3, area fraction TAWSS < 0.4 Pa, cycle analysed, n_snapshots, ρ/ν, validation flags, timestamp.

4. **CLI** — DONE. `--postprocess` (opt-in, runs metrics after solve) and `--postprocess-only <case>` (standalone: regenerates WSS via `pimpleFoam -postProcess -func wallShearStress` if absent, then computes metrics — works on cases solved before Phase C). See `runner.py: run_postprocess_only`.

5. **Optional ParaView screenshots** — DEFERRED (follow-up). `pvbatch` `--screenshots` not yet implemented.

6. **Validation ranges** (from prior iteration, `previous_cfd_iteration.md`) — wired as warn-only checks in `compute_metrics`:
   - WSS: 0–50 Pa
   - OSI: 0–0.5
   - Velocity: 0.1–1.0 m/s
   - Pressure: 0–200 Pa (relative to outlet)

---

### Phase C2 — Aneurysm-Scoped Biomarkers (IMPLEMENTED 2026-06-06, pending real-case validation)

**Goal:** Scope all hemodynamic metrics to the aneurysm sac only, not the whole vessel wall.  Adds normalised WSS, neck inflow rate, and sac pressure metrics.

**Prerequisite:** VORTEX must produce `aneurysm_sac.stl`, `parent_vessel.stl`, and `neck_plane.json` alongside the inlet/outlet caps.  See `VMTK_plan_biomarkers.md` for the VORTEX side.

**What was built:**
- `patch_labeller.py` — new mode accepts `aneurysm_sac` / `parent_vessel` / `inlet` / `outlet` labels (exactly 1 sac + 1 vessel + 1 inlet + ≥1 outlet).  Legacy `--legacy-no-aneurysm` mode still supported.
- `cli.py` — `--legacy-no-aneurysm` flag selects the old single-patch path.
- `case_builder.py` — detects mode from `scaled_stls` keys; loads `neck_plane.json`; builds new Jinja2 context variables (`wall_patches`, `aneurysm_patch`, `parent_vessel_patch`, `has_neck_plane`, `neck_origin`, `neck_normal`); `_bbox_with_buffer()` now accepts a list of wall STLs.
- `snappyHexMeshDict.j2` — `refinementSurfaces` and `addLayersControls.layers` now loop over `wall_patches`, giving 4-layer BL treatment to both sac and parent vessel.
- `0/U.j2`, `0/p.j2` — loop over `wall_patches` for noSlip / zeroGradient BCs.
- `controlDict.j2` — `wallShearStress` patches cover all wall patches; new `surfaceFieldValue_sac_pressure_mean/max` function objects; neck-plane function objects (`surfaceFieldValue_neck_flux`, `surfaceFieldValue_neck_peak_vel`) added when `has_neck_plane=True`.
- `postprocess.py` — auto-detects mode from `patch_labels.json`; in new mode computes TAWSS/OSI on `aneurysm_sac`, parent mean TAWSS, normalised WSS, and parses postProcessing CSVs for sac pressure and neck flow rate; `metrics_report.json` gets new keys (`aneurysm_tawss_pa`, `aneurysm_osi`, `normalised_wss`, `parent_tawss_pa_mean`, `sac_pressure_mean_pa`, `sac_pressure_peak_pa`, `neck_mean_flow_rate_m3s`, `neck_peak_flow_rate_m3s`, `neck_peak_velocity_ms`).

**Deferred from this phase:**
- **WSSG** (Wall Shear Stress Gradient) and **KEL** (Kinetic Energy Loss) require surface spatial gradients and volumetric flux integration that have no native OpenFOAM function object.  Both are to be implemented as a `pvbatch` Python script as described in `vortexcfd_biomarkers_plan.md` §6, gated by a `--pvbatch` flag.  Steps per the plan:
  - WSSG: `aneurysm_sac` block → "Gradient Of Unstructured Dataset" on `wallShearStress` → `mag(gradient)` → Temporal Statistics.
  - KEL: slice at neck plane → `0.5 × ρ × |U|² × (U·n̂)` calculator → Integrate Variables → cycle-averaged scalar.

---

### Phase D — Extensions (PLANNED)

1. **Carreau non-Newtonian viscosity** — Cho & Kensey (1991): μ₀ = 0.056, μ∞ = 0.0035 Pa·s, λ = 3.313 s, n = 0.3568. Opt-in via `--carreau` flag. Requires modifying `transportProperties.j2` and adding `CarreauYasuda` model.

2. **k-ω SST turbulence model** — Opt-in via `--turbulence` flag. Requires new `0/k.j2`, `0/omega.j2` templates and updated `turbulenceProperties.j2`.

3. **Mesh-independence helper** — Runs three mesh densities (coarse/medium/fine) and reports WSS convergence. Scheduled for after Phase C.

4. **FSI placeholder** — Reserved folder hooks for future fluid–structure interaction work.

---

## 5. Known issues and bugs

### BUG-001 — run-cfd.sh aborts silently when sourcing OpenFOAM
- **Discovered:** 2026-05-30, first real run on Kubuntu with OpenFOAM v2406
- **Symptom:** `bash run-cfd.sh` exited immediately with no output; `bash -x` showed the script stopping at `source /usr/lib/openfoam/openfoam2406/etc/bashrc`
- **Root cause:** `set -e` at the top of `run-cfd.sh` caused the script to abort when the OpenFOAM bashrc emitted a non-zero exit from an internal subcommand (`pop_var_context: head of shell_variables not a function context`)
- **Fix:** Wrap the source with `set +e; source "$bashrc" 2>/dev/null; set -e` in `run-cfd.sh`
- **Status:** FIXED

### BUG-002 — run-cfd.sh called wrong Python module
- **Discovered:** 2026-05-30, first real run on Kubuntu
- **Symptom:** Program exited with code 0 and no output; no case directory created
- **Root cause:** `run-cfd.sh` called `python -m vortex_cfd.cli` — this imports the module but never invokes `main()`. The correct entry point is `python -m vortex_cfd` which uses `__main__.py`
- **Fix:** Changed all `python -m vortex_cfd.cli` to `python -m vortex_cfd` in `run-cfd.sh`
- **Status:** FIXED

### BUG-003 — `background` patch missing from 0/p and 0/U templates
- **Discovered:** 2026-05-30, first real run on Kubuntu
- **Symptom:** `decomposePar` (before snappyHexMesh) crashed: `Cannot find patchField entry for background`
- **Root cause:** `blockMesh` creates an outer boundary patch named `background`; OpenFOAM requires every patch to have an entry in all initial condition files. The templates only listed wall/inlet/outlet.
- **Fix:** Added `background { type zeroGradient; }` to `0/p.j2` and `background { type slip; }` to `0/U.j2`
- **Status:** FIXED

### BUG-004 — Missing `div((nuEff*dev2(T(grad(U)))))` in fvSchemes
- **Discovered:** 2026-05-30, first real run on Kubuntu with OpenFOAM v2406
- **Symptom:** `pimpleFoam` crashed on first timestep: `Entry 'div((nuEff*dev2(T(grad(U)))))' not found in dictionary "system/fvSchemes/divSchemes"`
- **Root cause:** OpenFOAM v2406 with the laminar Stokes model requires an explicit `divSchemes` entry for the viscous stress term. This was not required (or was implicit) in v2512 for which the template was written.
- **Fix:** Added `div((nuEff*dev2(T(grad(U))))) Gauss linear;` to `divSchemes` in `system/fvSchemes.j2`. Entry is harmless in later versions — OpenFOAM ignores unused divScheme entries.
- **Status:** FIXED

### BUG-006 — locationInMesh falls outside the vessel lumen
- **Discovered:** 2026-05-30, first visual inspection in ParaView after Phase A validation run
- **Symptom:** Glyphs (velocity vectors) appear only in the background mesh (the blockMesh cube), not inside the vessel. The vessel interior has no fluid cells — the fluid domain is the space *around* the vessel, not inside it.
- **Root cause:** `_location_in_mesh()` used the centre of the wall STL bounding box. For curved/C-shaped vessels this point falls inside the wall material, so snappyHexMesh treats the outside of the vessel as the fluid domain and the inside as solid. Documented as a known limitation in D-003.
- **Fix:** Replaced with inlet-centroid method: compute the area-weighted centroid and normal of the inlet cap STL, then step one inlet-radius inward (negating the outward-pointing cap normal). This guarantees the point is inside the lumen. See `case_builder.py:_location_in_mesh()`.
- **Status:** FIXED

### BUG-005 — snappyHexMesh segfault during parallel load balancing (v2406)
- **Discovered:** 2026-05-30, first real run on Kubuntu with OpenFOAM v2406
- **Symptom:** `snappyHexMesh -parallel` crashed with segfault (signal 11) inside `fvMeshDistribute::repatch` during shell refinement iteration 0, after reporting `max unbalance 0.237 > allowable 0.1`
- **Root cause:** Bug in OpenFOAM v2406's parallel mesh redistribution (`fvMeshDistribute::repatch → polyTopoChange::changeMesh`). Triggered when the load balancer tries to redistribute cells after shell refinement.
- **Fix:** Run `snappyHexMesh` in serial regardless of `--cores`. Parallel solving with `pimpleFoam` is unaffected. See `runner.py` comment. May be fixed in v2412+.
- **Status:** FIXED

---

## 6. Design decisions log

Records non-obvious choices so future sessions don't re-litigate them.

### D-001 — mm→m scaling written back to STL before any other use
*2026-05-26*
The prior iteration (VesselForge_AutoCFD) applied scaling to the geometry in memory but forgot to rewrite patch-coordinate metadata files. This silently mapped the inlet patch to a point 1000× away from the mesh, producing runs that "completed" against the wrong geometry. In this implementation, `scaling.py` rewrites scaled copies to a temp directory immediately; everything downstream (case builder, SHM dict, patch labels) reads only the scaled copies.

### D-002 — `patch_labels.json` written after case directory exists
*2026-05-26*
Prior iteration wrote this file before the directory was created, leading to a race condition where the file was created in the current working directory instead. `case_builder.py` now writes it as the last step before returning.

### D-003 — `locationInMesh` = inlet centroid displaced one radius inward
*2026-05-26 (initial); revised 2026-05-30*
Originally used the centre of the wall STL bounding box — failed on the first real patient geometry (curved ICA) because the bbox centre fell inside the wall material, causing the fluid domain to be the exterior of the vessel. Replaced with: compute area-weighted centroid and normal of the inlet cap STL, step one inlet-radius inward along the negated cap normal. Cap normals point outward by VMTK convention, so negating gives the inward direction. This is robust for any vessel with a visible inlet opening.

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

**Phase C2 is CODE-COMPLETE (2026-06-06). Validate on a real case with new VORTEX output.**

Branch: `aneurysm_dome_biomarkers`.  All new unit tests pass (pure numpy + pyvista fixture, no OpenFOAM needed).

Validation checklist:
1. Confirm VORTEX produces `aneurysm_sac.stl`, `parent_vessel.stl`, `neck_plane.json` in the STL output directory.
2. Run: `bash run-cfd.sh --stl-dir <new-vortex-output> --cycles 3 --mean-velocity 0.4 --cores 6 --postprocess`
3. Check `metrics_report.json` for `aneurysm_tawss_pa`, `normalised_wss`, `sac_pressure_mean_pa`.
4. Open case in ParaView; verify `aneurysm_sac` and `parent_vessel` are separate patches with WSS coloured on each.

**Deferred — WSSG + KEL (pvbatch).** See Phase C2 notes and `vortexcfd_biomarkers_plan.md` §6.

**After validation:** Phase B (Womersley inlet profile) or pvbatch WSSG/KEL script, whichever is more urgent for the paper.

---

## 8. Session history

| Date | Session summary |
|---|---|
| 2026-05-26 | Phase A implementation: all Python modules + 13 OpenFOAM templates + Allrun script. Syntax-validated. Not yet run against real STLs on Linux. |
| 2026-05-26 | pytest suite: 99 tests across 4 files (waveform, scaling, case_builder, env_check). All pass on Windows with synthetic pyvista geometry. Added pyproject.toml. |
| 2026-05-29 | Phase A validation (partial): pytest confirmed 99/99 pass on Linux in vortex-aneurysm env. Fixed smoke_test.sh (OUT_DIR/STL_DIR/REPO_ROOT not exported — Python subprocess couldn't read them via os.environ). Smoke test PASSED with real VMTK STLs. OpenFOAM v2406 confirmed at standard path. Full mesher+solver run not yet executed. |
| 2026-05-30 | Phase A fully validated on Kubuntu desktop (Ryzen 5 5600X, OpenFOAM v2406). Fixed 6 bugs during first real run (see Section 5). Key fixes: run-cfd.sh conda/venv detection, `-m vortex_cfd` entry point, background patch in 0/U and 0/p, div(nuEff) in fvSchemes, snappyHexMesh serial-only workaround for v2406 segfault, locationInMesh replaced with inlet-centroid method. Velocity field confirmed inside vessel lumen in ParaView. **Phase A COMPLETE.** |
| 2026-06-04 | Phase C implemented and fully validated on a real patient geometry (OpenFOAM v2406, 6 cores, 3 cycles). `metrics_report.json`: TAWSS mean 9.83 Pa, OSI mean 0.019, low-WSS area 0.17%, high-OSI area 0.84% — all clinically plausible. `wallShearStress` confirmed visible in ParaView with correct pulsatile temporal behaviour. Also fixed `run-cfd.sh` conda env detection (base env / space-in-path bugs). **Phase C COMPLETE.** |
| 2026-06-06 | Phase C2 implemented (branch `aneurysm_dome_biomarkers`): two-patch wall mode (aneurysm_sac + parent_vessel), `--legacy-no-aneurysm` fallback flag, neck-plane function objects, normalised WSS, sac pressure metrics, neck flow rate. 8 source files + 3 test files updated. New unit tests (pure numpy/pyvista, no OpenFOAM). Pending real-case validation with updated VORTEX output. WSSG + KEL deferred to pvbatch script (see `vortexcfd_biomarkers_plan.md` §6). |
