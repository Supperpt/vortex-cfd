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
| Inlet waveform | Ford et al. (2005) ICA archetype | Literature-validated default, normalised to mean=1; scaled by `--mean-velocity`. Bundled `data/ica_ford2005.csv`. Override via `--waveform`. |
| Inlet BC | flowRateInletVelocity (default), Womersley via `--womersley` | Parabolic by default; the exact Womersley profile is opt-in (Phase C4, D-007) |
| Outlet BC | inletOutlet (U), fixedValue 0 (p) | Prevents recirculation instability |
| Wall BC | noSlip (U), zeroGradient (p) | Standard rigid-wall |
| BL layers | 4 prismatic, expansion 1.3, finalLayerThickness 0.3 | Required for WSS accuracy |
| BL refinement | level (3 4) → ~0.125 mm at wall | From 2 mm background cells |
| checkMesh abort | maxNonOrtho > 70°, maxSkewness > 20 | Cells that crash pimpleFoam later. Skewness 20 = OpenFOAM `maxBoundarySkewness` (internal-4 is the meshing target). See BUG-009. |

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
│   ├── waveform.py                ← loads default Ford ICA waveform (data/) + user CSV
│   ├── data/
│   │   ├── ica_ford2005.csv       ← default inlet waveform (Ford et al. 2005, mean=1)
│   │   └── generate_ica_ford2005.py  ← reproducible generator (digitised Table 2)
│   ├── case_builder.py            ← geometry analysis + Jinja2 rendering
│   ├── womersley.py               ← Phase C4: analytical Womersley inlet + boundaryData writer
│   ├── neck.py                    ← Phase C3: neck orifice fitting + inflow metrics
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
- `scaling.py`: detects mm coordinates (bounding-box *extent* > 1.0, translation-invariant — see BUG-014) and scales ×0.001 before saving to temp dir; writes canonical names (wall.stl, inlet.stl, outlet_0.stl, …)
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

### Phase C2 — Aneurysm-Scoped Biomarkers (COMPLETE 2026-06-13)

**Validated on two real patients** (AA_010, AA_004; 0.30 m/s ICA inlet, 3 cycles, 6 cores). The WSS/OSI biomarkers are trustworthy and physiologically plausible:
- Aneurysm TAWSS mean 3.1 / 3.8 Pa, parent TAWSS mean 6.0 / 6.9 Pa, normalised WSS 0.52 / 0.56 (sac < parent, as expected).
- Aneurysm OSI mean 0.023 / 0.008 (max < 0.39); low-WSS area 0% / 8%. Validation flags (TAWSS in 0–50, OSI in 0–0.5) pass. Last cycle analysed, 50 snapshots, full 3 cycles run.

**The neck inflow / peak-velocity metrics are NOT trustworthy yet** — three issues found during validation, deferred to **Phase C3** (see below). Sac pressure (≈4–14 mmHg gauge) is plausible but exceeds the old 0–200 Pa note (that note is too conservative; pressure is not in the enforced validation flags).

(Original implementation notes below, IMPLEMENTED 2026-06-06.)

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
- **WSSG** (Wall Shear Stress Gradient) and **KEL** (Kinetic Energy Loss) require surface spatial gradients and volumetric flux integration that have no native OpenFOAM function object.  Both are to be implemented as a `pvbatch` Python script as described in `docs/planning/vortexcfd_biomarkers_plan.md` §6, gated by a `--pvbatch` flag.  Steps per the plan:
  - WSSG: `aneurysm_sac` block → "Gradient Of Unstructured Dataset" on `wallShearStress` → `mag(gradient)` → Temporal Statistics.
  - KEL: slice at neck plane → `0.5 × ρ × |U|² × (U·n̂)` calculator → Integrate Variables → cycle-averaged scalar.

---

### Phase C3 — Neck inflow metrics (IMPLEMENTED 2026-08-04, awaiting validation)

Found during C2 validation of AA_010/AA_004; BUG-010, BUG-011 and CAVEAT-012 all
traced back to the same root cause, so all three are closed by one redesign.

**The decision: the `surfaceFieldValue` function objects were deleted, not fixed.**
The clinical quantity is the *neck inflow rate* — the volume entering the sac per
second, i.e. the integral of only the **inward** part of U·n̂ over the orifice. No
`surfaceFieldValue` operation can express a positive-part integral;
`areaNormalIntegrate` gives the *net* flux, which for a sealed sac averages to ~0
over a cycle by conservation. Chasing the FO fixes would have produced a correct
implementation of the wrong quantity. The metrics are now computed in Python
(`vortex_cfd/neck.py`) by slicing the saved `U` volume snapshots at the orifice.

**The orifice comes entirely from `aneurysm_sac.stl`.** VORTEX clips the sac at the
neck, so its open boundary loop *is* the orifice; an SVD fit to that loop gives
origin, normal and extent that are mutually consistent by construction (measured
planarity ratio 4.8e-7 on a clipped sphere). Notably the normal's **sign is derived
geometrically**, by orienting it towards the sac's area-weighted centroid — *not*
taken from `neck_plane.json`, whose sign convention is undocumented and whose normal
`case_builder` never normalised. A flipped normal would have reported the **outflow**
rate: plausible magnitude, right units, wrong metric, and undetectable without
ParaView. That is the same failure class as BUG-011. `neck_plane.json` is kept purely
as a cross-check and is recorded in the output.

Radius is `r_eff = sqrt(A_loop/π)` rather than `r_max`, deliberately biasing the disc
*small*: over-inclusion re-creates CAVEAT-012 (the disc reaches into the parent
vessel), while under-inclusion is only a mild area bias.

**Status: opt-in and unvalidated.** No solved case was available when this landed, so
`--neck-metrics` defaults off and `metrics_report.json` still emits
`"neck_metrics": "disabled_pending_validation"` by default — the released output is
byte-identical to v1.0.0. See the flip-the-switch procedure below.

**Self-validation without ParaView.** The report includes `net_to_inflow_ratio`. For a
sealed sac the cycle-mean *net* flux must be ~0 (physics, not convention), so a ratio
near 1 proves the disc is still cutting the parent vessel. CAVEAT-012 is therefore now
self-detecting rather than needing to be spotted by eye.

#### Flip-the-switch procedure (after validating on a real case)

1. `--postprocess-only <case> --neck-metrics` on any solved aneurysm case.
2. Check `validation.net_flux_near_zero` in `metrics_report.json`. If false, the disc
   is over-reaching.
3. In ParaView: `Slice` on `internalMesh` at the `origin`/`normal` from
   `neck_plane_resolved.json`, then `Clip → Sphere` at that file's `radius_m`. Confirm
   the disc covers the sac orifice and **no parent-vessel lumen**. If it over-reaches,
   hand-edit `radius_m` in `neck_plane_resolved.json` and re-run `--postprocess-only`
   — no re-solve needed.
4. Flip `# PHASE-C3-VALIDATION-SWITCH` in `cli.py` from `default=False` to
   `default=True`, change `"status": "experimental_unvalidated"` to `"validated"` in
   `postprocess._neck_report_block`, and drop the EXPERIMENTAL notes from README.

(Then resume the **WSSG + KEL pvbatch** work deferred from C2 — see
`docs/planning/vortexcfd_biomarkers_plan.md` §6. KEL is now one line over the same
clipped disc: `0.5·ρ·|U|²·(U·n̂)`.)

---

### Phase C4 — Womersley inlet profile (IMPLEMENTED 2026-08-04, awaiting validation)

Opt-in `--womersley` flag replacing the default parabolic `flowRateInletVelocity`
inlet with the exact Womersley profile: FFT of the waveform, then a complex-Bessel
radial shape per harmonic, each **area-normalised** so the area-weighted mean of
U(r,t) reproduces `U_mean · waveform(t)`. That identity is what the unit tests assert
(verified to 1.0000 by radial quadrature, ±2 % across the cycle from truncating at 8
harmonics). Delivered via `timeVaryingMappedFixedValue` boundaryData — see D-007 for
why not `codedFixedValue`.

Re-applied from the abandoned `phase_b_womersley_inlet` branch (forked pre-Phase-C2,
~18 commits stale — merging it would have regressed shipped work). Deliberately
**dropped** what that branch bundled: the `print`→`logging` migration (would have
deleted the `run_log.md` feature and the BUG-009 checkMesh regex fallback),
`--dry-run`, physiological flow-rate validation, and the snappyHexMesh retry. Design
notes in `docs/planning/fable_womersley_plan.md`.

**The prototype's real bug, now fixed.** `timeVaryingMappedFixedValue` has no
equivalent of `flowRateInletVelocity`'s `outOfBounds repeat`: past the last supplied
time it *holds the final value* rather than wrapping. The branch wrote a single
cycle, so any `--cycles > 1` run would have frozen the inlet at the cycle-1 end value
for every subsequent cycle — a silent, physically plausible-looking failure.
`womersley_velocities` now takes `cycles` and tiles the normalised cycle across the
whole run, plus a terminal sample at exactly `cycles·T`. `TestMultiCycleCoverage`
pins this.

**Implementation notes.**
- `_location_in_mesh` was promoted to `_inlet_geometry`, which returns the centroid,
  inward normal, equivalent radius and interior point it already computed internally.
  `_inlet_area` is retained for the parabolic flow table — the equivalent radius is
  only for the Womersley radial shape.
- `build_case` now returns `(case_dir, inlet_params)`.
- boundaryData is written between `snappyHexMesh` and `checkMesh` (inlet face centres
  only exist after meshing; files must precede `decomposePar`). Unlike the neck
  metrics, this step **aborts on failure** — a missing boundaryData directory leaves
  the solver with an unreadable inlet BC, so failing early costs seconds rather than
  hours.
- Near-wall retrograde flow is deliberately not clipped; it is physically correct for
  Womersley flow, and clipping would inflate the area-mean exactly at diastole.
- Cost: a 3-cycle run over a 1200-face inlet writes 301 time directories, ~7 MB.

**Not yet run against a real OpenFOAM case.** Validation step 2 in NEXT ACTION.

---

### Phase D — Extensions (PLANNED)

1. **Carreau non-Newtonian viscosity** — Cho & Kensey (1991): μ₀ = 0.056, μ∞ = 0.0035 Pa·s, λ = 3.313 s, n = 0.3568. Opt-in via `--carreau` flag. Requires modifying `transportProperties.j2` and adding `CarreauYasuda` model.

2. **k-ω SST turbulence model** — Opt-in via `--turbulence` flag. Requires new `0/k.j2`, `0/omega.j2` templates and updated `turbulenceProperties.j2`.

3. **Mesh-independence helper** — Runs three mesh densities (coarse/medium/fine) and reports WSS convergence. Scheduled for after Phase C.

4. **FSI placeholder** — Reserved folder hooks for future fluid–structure interaction work.

---

### Phase E — Extended rupture-risk biomarker suite (PLANNED, opened 2026-08-03)

Goal: get the pipeline to a place where it can compute the full set of biomarkers
identified by literature review as likely discriminants for **small (<5mm) intracranial
aneurysm rupture**, beyond the current TAWSS/OSI/normalised-WSS/sac-pressure set. See
**`docs/planning/Biomarcadores_candidatos.md`** for the full research summary
(evidence, effect sizes, and an overfitting caution for small-aneurysm cohorts).

Candidate metrics, roughly in order of expected significance for this population:
- **Size Ratio (SR)** — morphological, aneurysm/parent-vessel dimension ratio.
  **Out of scope for vortex-cfd** — pure STL/mesh geometry metric, no CFD dependency.
  Tracked upstream: `Supperpt/VORTEX` milestone "Compute morphological rupture-risk
  biomarkers (SR, NSI, UI)" (#2).
- **OSI** — already implemented (Phase C2); flagged in the research as the primary
  size-specific discriminant, so existing validation should hold up.
- **High Shear Concentration Ratio (HSCR)** and **Wall Shear Stress Divergence
  (WSSD)** — local WSS concentration/traction metrics, computable from the existing
  `wallShearStress` field on `aneurysm_sac` (likely a `pvbatch` gradient job, same
  family as the deferred WSSG work in Phase C2/C3).
- **Oscillatory Velocity Index (OVI)** — 3D flow-instability analogue of OSI, needs a
  velocity-field (not just wall-field) time-series analysis.
- **Non-Sphericity Index (NSI) / Undulation Index (UI)** — morphological. **Out of
  scope for vortex-cfd**, same as SR — tracked in `Supperpt/VORTEX` milestone #2.
- **Flow Complexity Ratio (FCR)** — qualitative jet-concentration/vortex-count metric;
  needs a concrete quantitative definition before implementation.

**Scoping decided (2026-08-03):** SR, NSI, and UI are pure mesh/STL geometry metrics
with no CFD dependency, so they're computed upstream in `Supperpt/VORTEX` (milestone
#2) rather than in vortex-cfd, to avoid duplicating geometry logic across the two
projects. This phase's scope is therefore HSCR, WSSD, OVI, and FCR.

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

### BUG-007 — Phase C2 `surfaceFieldValue` function objects use invalid ESI v2406 syntax
- **Discovered:** 2026-06-10, first real run of the Phase C2 two-patch pipeline (OpenFOAM v2406, 6 cores)
- **Symptom:** `pimpleFoam` aborted on startup: `Entry 'writeFields' not found in dictionary ".../surfaceFieldValue_sac_pressure_mean"`
- **Root cause:** The `surfaceFieldValue` blocks in `controlDict.j2` were written with keys that ESI v2406 does not accept: `surfaceType`/`patches` (correct: `regionType`/`name`), no `writeFields` entry (mandatory), and a non-existent `surfaceType plane; basePoint; normalVector` form for the neck plane. These blocks are new in Phase C2 and had never reached the solver before — Phase C (2026-06-04) only emitted `wallShearStress` + `fieldAverage`, and the Phase C2 unit tests only assert that the template renders the FO *names*, not that OpenFOAM accepts the dict syntax.
- **Fix:** Rewrote the FOs against the installed v2406 caseDicts/source:
  - Patch FOs: `regionType patch; name <patch>; writeFields false;`
  - Neck FOs: `regionType sampledSurface; name ...; sampledSurfaceDict { type plane; planeType pointAndNormal; pointAndNormalDict { point; normal; } source cells; interpolate true; }`
  - Neck flux operation `sum` (meaningless on U) → `areaNormalIntegrate` (true volumetric flow rate ∫U·n dA).
- **Status:** FIXED

### BUG-008 — Phase C2 post-process parser read the wrong `.dat` filename
- **Discovered:** 2026-06-10, alongside BUG-007 (found by reading the v2406 `writeFile`/`fieldValue` source)
- **Symptom:** Would have silently returned no rows (sac pressure / neck metrics → `None`) even on a successful solve.
- **Root cause:** `postprocess.py:_read_surface_field_value` read `surface_fieldValue.dat`; ESI names the file after the FO `typeName`, i.e. `surfaceFieldValue.dat`.
- **Fix:** Corrected the filename in `postprocess.py` and the test fixture in `tests/test_postprocess.py`.
- **Status:** FIXED

### BUG-009 — checkMesh skewness gate too strict (false abort on good mesh)
- **Discovered:** 2026-06-13, batch of 3 patient runs (AA_010, AA_004, AA_011). AA_011 was run #3.
- **Symptom:** AA_011 meshed cleanly (snappyHexMesh "Finished meshing without any errors", layers added on both wall patches) but the run aborted at the post-mesh gate: `ERROR: maxSkewness = 6.66 > 4 threshold`. Solver never ran; no `metrics_report.json`. AA_010/AA_004 (runs #1/#2) had already completed (chained with `&&`).
- **Root cause:** `runner._check_mesh_quality` hard-aborted on checkMesh's global "Max skewness" > 4. That is OpenFOAM's *internal*-face meshing target, not a fatal limit. snappyHexMesh's own quality check passed with 0 faces over its internal-4/boundary-20 criteria (log line 868); checkMesh re-flags boundary-layer faces against a single threshold of 4. AA_011 had 29 skewed faces / 3.12 M (0.0009%), max 6.66, non-ortho 69.96 (passed), 0 negative volumes — a usable mesh.
- **Fix:** Raised the skewness abort threshold 4 → 20 in `runner._check_mesh_quality`, matching OpenFOAM `maxBoundarySkewness`. Non-orthogonality (70°) unchanged; negative volumes/determinant remain the genuinely fatal checks. Skewness is a cell-*shape* metric (independent of cell size — not a "too fine" issue); the deeper lever for a genuinely skewed case is STL smoothing/remeshing in VORTEX. Literature: OpenFOAM snappyHexMesh guide + CFD Support (skewness ≤ 20 usable); aneurysm-CFD robustness study (Steinman group, PMC5469453) — bulk WSS robust, OSI/WSSG resolution-sensitive, so keep the skewed-face count small.
- **Also fixed (found while patching the gate):** the non-orthogonality regex `r"Max non-orthogonality\s*=\s*..."` never matched checkMesh's actual output (`Mesh non-orthogonality Max: 69.96 average: ...`), so the 70° gate had *silently never fired* — a mesh with non-ortho 85 would have passed. Corrected to `r"non-orthogonality Max:\s*([\d.]+)"`. Regression-tested in `tests/test_runner.py`.
- **Status:** FIXED

### BUG-010 — `neck_peak_velocity_ms` always null (vector parser)
- **Discovered:** 2026-06-13, C2 validation of AA_010/AA_004 (both report `null`).
- **Symptom:** `neck_peak_velocity_ms` is `null` in every report.
- **Root cause:** the `neck_peak_vel` FO writes a parenthesised vector per line, `<t>\t(vx vy vz)`. `postprocess._read_surface_field_value` does `line.split()` → `parts = ['<t>', '(vx', 'vy', 'vz)']`, then `float('(vx')` raises `ValueError`, which the `except ValueError: pass` swallows → every row dropped → empty series → `None`. Confirmed by reproduction. Secondary: the FO operation `max(U)` is component-wise (a vector), not the peak speed.
- **Fix (2026-08-04, C3):** parser strips `()` before parsing and now **counts malformed rows and warns**. The silence was the more serious half of this bug: a total parse failure surfaced only as an unexplained `null`. The FO itself is gone — peak velocity comes from the Python disc sampler. A lesson worth keeping: the old test `test_vector_field_returns_magnitude` fed `3.0 4.0 0.0`, a format OpenFOAM never writes, so it passed green while the real path was 100 % broken.
- **Status:** FIXED

### BUG-011 — `neck_peak_flow_rate_m3s` uses sign-naive max
- **Discovered:** 2026-06-13, C2 validation. AA_010 reports peak (−1.33e-6) smaller in magnitude than the mean (−1.95e-6) — impossible for a true peak.
- **Root cause:** `postprocess.py` computed `float(max(vals))` on the *signed* flux series. For an outward-pointing neck normal the flux is negative all cycle, so `max()` returned the least-negative (smallest-magnitude) value. AA_004 was correct only by luck (positive flux).
- **Fix (2026-08-04, C3):** obsolete by redesign. The reported quantity is now a positive-part integral, `Σ A·max(U·n̂, 0)`, which cannot return a spurious peak from an all-negative series; `tests/test_neck.py::test_all_outflow_gives_zero_inflow` locks that in. The underlying hazard — a sign convention that silently inverts a metric — is addressed structurally by deriving the normal's sign from the sac geometry instead of from `neck_plane.json`.
- **Status:** FIXED

### CAVEAT-012 — neck flux integrates the whole vessel cross-section
- **Discovered:** 2026-06-13, C2 validation. Neck mean flux (~1.9–2.5e-6 m³/s) ≈ ICA parent throughput (v·A ≈ 0.30 × ~9e-6 ≈ 2.7e-6), not a small sac-neck inflow (which should net ~0 over a cycle).
- **Root cause:** the neck FO used an *infinite* `sampledSurface` plane; `areaNormalIntegrate(U)` summed over the entire plane∩fluid intersection, which extended outward cuts the parent vessel lumen — so it measured vessel throughput, not sac-neck inflow.
- **Fix (2026-08-04, C3):** the sampling surface is clipped to a disc of radius `r_eff` about the fitted orifice centroid, using `clip_scalar` on an in-plane radius (measured −1.3 % area error, vs +3.4 % for a staircased cell-centre mask). The evidence that first revealed this — neck flux ≈ parent throughput — is now an automated check: `validation.net_to_inflow_ratio` near 1 means the disc still over-reaches.
- **Status:** FIXED IN CODE, awaiting ParaView confirmation on a real case (no solved case available at implementation time).

### BUG-013 — pimpleFoam GAMG pressure-solver FPE on rough (MRI-derived) meshes — DEFERRED (Phase B)
- **Discovered:** 2026-06-14, re-run of AA_011 (the only MRI-derived case; larger voxel → rougher STL).
- **Symptom:** `pimpleFoam` aborted with a floating-point exception in the GAMG pressure solver (`GAMGSolver::scale`) on the **first** PIMPLE iteration at Courant 0 (exit 136). Velocity solved fine; only the pressure linear solve failed. Captured in the case `run_log.md`.
- **Root cause:** not an impossible case — checkMesh passed (0 negative volumes, 0 illegal faces, non-ortho 70, determinant OK), but boundary-layer addition was distressed (only 74% layer coverage on the sac), leaving sliver cells. GAMG agglomeration produces a singular coarse-grid coefficient on slivers → divide-by-zero → FPE.
- **Proposed answer (not implemented):** pressure solver GAMG → PCG/DIC in `fvSolution.j2` (identical converged WSS, ~1.3–1.5× slower, no agglomeration step to go singular); best as an **auto-retry fallback** in `runner.py` so GAMG stays the fast default and unattended batches self-heal.
- **Decision:** DEFERRED — we are in a proof-of-concept sampling phase, and the rough-mesh root cause is being handled upstream in VORTEX remesh (milestone "Conclude remesh functionality", VORTEX#3/#4/#5). Tracked as `Supperpt/vortex-cfd#5`, milestone Phase B. Revisit if remesh alone proves insufficient.
- **Status:** OPEN (Phase B, deferred)

### BUG-014 — mm/m detection used distance-from-origin, not extent (translation-variant)
- **Discovered:** 2026-06-23, code review.
- **Symptom:** none observed on real cases (see below) — flagged as a latent correctness risk.
- **Root cause:** `scaling._any_in_mm` tested `max(abs(v) for v in mesh.bounds) > 1.0`, i.e. the largest *absolute coordinate* (distance from the origin), not the mesh's *extent*. DICOM coordinates are tied to the scanner isocenter, so a metre-scale mesh offset from (0,0,0) could in principle be misjudged. **Severity reassessed: NOT critical for this project.** For cerebral-aneurysm head scans the anatomy sits inside the ~0.35 m scanner bore near isocenter, so metre coordinates never approach 1.0 m and mm coordinates (vessels ~5–50 mm) are always ≫ 1.0 — the origin-distance and extent checks return the *identical* verdict on every real VORTEX case. The review's "catastrophic, metre mesh at Z=1.5 m" scenario cannot arise for head DICOM. Fixed as defensive hardening, not a live-bug fix. (Note: `_any_in_mm` was unchanged since the first commit — it was never "altered"; that memory was the `locationInMesh` bbox→inlet-centroid change, BUG-006/D-003.)
- **Fix:** Compute `max_extent = max(xmax-xmin, ymax-ymin, zmax-zmin)` and compare to `_BBOX_THRESHOLD`, making detection translation-invariant. Regression test `test_metre_geometry_offset_from_origin_not_detected_as_mm` (metre sphere offset +1.5 m → not flagged) added in `tests/test_scaling.py`.
- **Status:** FIXED

### BUG-015 — `--postprocess-only` silently assumed 3 cycles (wrong analysis window)
- **Discovered:** 2026-07-02, code review (Fable).
- **Symptom:** `--postprocess-only <case>` always analysed the window `t >= 2·T` regardless of how many cycles the case was actually solved with. A case solved with 5 cycles got cycles 3–5 averaged (not just the last); a 2-cycle case (`endTime = 2·T`) got only the final snapshot, giving degenerate stats (OSI ≈ 0 everywhere).
- **Root cause:** `compute_metrics` already supports `cycles=None` to derive the last cycle from the case's own snapshot times (`t_start = t_max − T`), and its docstring documents that path for standalone post-processing. But it was **unreachable from the CLI**: `--cycles` had `default=3`, so Click always supplied `3`, and the post-process-only branch forwarded it unconditionally (`cli.py: run_postprocess_only(..., cycles=cycles)`). The auto-derive branch was dead code.
- **Fix:** `--cycles` now defaults to `None`. A full run substitutes `3` (`if cycles is None: cycles = 3`) after the post-process-only branch, so `--postprocess-only` without `--cycles` forwards `None` and `compute_metrics` derives the window from the data; passing `--cycles N` still pins it. Regression tests in `tests/test_cli.py` (new file): post-process-only forwards `None` when omitted / the explicit value when given, and a full run still defaults to 3.
- **Status:** FIXED

### BUG-016 — `run-cfd.sh` hardcoded `~/miniconda3` for the conda-activation fallback
- **Discovered:** 2026-07-02, code review (Fable).
- **Symptom:** In the fallback branch (active env is not `vortex-aneurysm` and no `.venv`), `run-cfd.sh` sourced `$HOME/miniconda3/etc/profile.d/conda.sh`. On this machine conda lives at `~/miniforge3`, so the launcher aborted with "Could not find conda" even though conda was installed and working.
- **Fix:** Discover the conda base instead of hardcoding a distribution: prefer `conda info --base` when `conda` is on PATH, then fall back to a list of common locations (`~/miniforge3`, `~/mambaforge`, `~/miniconda3`, `~/anaconda3`, `/opt/conda`). `bash -n` syntax-checked.
- **Status:** FIXED

### Cleanup — postprocess polish (2026-07-02, code review Fable)
Three minor `postprocess.py` items from the same review, fixed as hardening (no observed live failure):
- **Dead validation branch.** `in_wss_range` tested `stats["tawss_pa"]["min"] >= 0.0`, always true (TAWSS is a mean of vector magnitudes → non-negative by construction). Reduced to the meaningful upper-bound check only.
- **Coupled parent-vessel weights.** In new mode the parent TAWSS reused the sac's `times` for its trapezoidal weights (`w_pv = _time_weights(times)`), discarding the parent read's own returned times. Now uses `times_pv` from the parent read — identical for a normal solve, but no longer coupled to the sac read.
- **`_read_surface_field_value` only read the earliest time dir.** A restarted solve writes one time-named subdir per restart; the parser took `time_dirs[0]` and ignored later ones, so a restart would drop the last-cycle rows. Now merges all time subdirs in numerical order (later dir wins on overlapping boundary times). Currently only exercises sac pressure (neck metrics disabled), but matters when those return. Regression test `test_merges_multiple_restart_time_dirs` in `tests/test_postprocess.py`.

### Feature — `run_log.md` per-case run log (2026-06-13)
`runner.py` now writes a markdown `run_log.md` into each case directory: a header (case path, start time, pipeline) plus one section per pipeline step (label, command, full stdout/stderr in a fenced block) and the mesh-quality verdict, written **incrementally** so the log survives a mid-run `sys.exit` abort. Motivated by BUG-009 diagnosis being hampered by no persisted solver log (output was stdout-only). Helpers `_log_init`/`_log_append`; no-ops safely if the file can't be written. Tested in `tests/test_runner.py`.

### Note — Python 3.9 compatibility (2026-06-10; superseded 2026-07-02)
Originally the shared `vortex-aneurysm` conda env was pinned to Python 3.9 (vtk/VMTK dependency constraints blocked an in-place upgrade). The code used 3.10+ `X | Y` union type hints, so the env could not import it (`run-cfd.sh` failed at `import click`). `from __future__ import annotations` was added to all `vortex_cfd/*.py` modules and `requires-python` set to `>=3.9`.
**Superseded:** the env is now Python 3.10 (`vortex-aneurysm` reports 3.10.16) and `pyproject.toml` sets `requires-python = ">=3.10"` with 3.10–3.12 classifiers. The `from __future__ import annotations` lines are kept (harmless, and cheap insurance against a future 3.9 env). vortex-cfd does **not** depend on VMTK — it only reads VORTEX's STL output — so the two tools can live in separate envs if a future Python bump is wanted.

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

### D-006 — Default STL naming scheme + auto-labelling (unattended batch runs)
*2026-06-13*
To run several cases unattended, `patch_labeller.py` now auto-labels from filenames before falling back to the interactive prompt. Scheme (case-insensitive stems): `aneurysm.stl`→aneurysm_sac, `wall.stl`→parent_vessel, `inlet.stl`→inlet, `outlet_<N>.stl`→outlet (legacy mode: `wall.stl`→wall). Auto-labelling fires **only if every STL in the dir matches the scheme and validation passes**; any unrecognised filename or count mismatch falls back to prompting for *all* files (user's choice — no partial/heuristic guessing). Mode detection is unaffected: `scaling.scale_stls` keys the canonical dict off the *label*, so `wall.stl`→`parent_vessel` key in new mode, and `case_builder`'s `is_legacy = "wall" in scaled_stls` only trips in true legacy mode. `neck_plane.json` is auto-loaded separately (not labelled) and optional. See `patch_labeller._auto_label` / `_label_from_filename`; tests in `tests/test_patch_labeller.py`.

### D-005 — `pimpleFoam` with `nOuterCorrectors 2`
*2026-05-26*
The PIMPLE algorithm with 2 outer correctors gives a good balance between stability and cost at Co < 1. Increasing to 3 would be safer for very coarse meshes but adds 50% cost per timestep.

### D-007 — Womersley via `timeVaryingMappedFixedValue`, not `codedFixedValue`
*2026-08-04*
The original README/LLM.md wording implied a `codedFixedValue` inlet with the Fourier–Bessel series evaluated inside an OpenFOAM coded BC. It avoids writing boundaryData files and covers all simulation time automatically. **Rejected**: it needs runtime C++ compilation (`dynamicCode`, fragile across ESI versions and often unavailable on locked-down HPC nodes); OpenFOAM has no native complex Bessel function, so `J0` of a complex argument would have to be hand-coded in C++; and the maths would stop being testable in Python. `timeVaryingMappedFixedValue` needs no compilation, keeps the maths in `scipy`, and is unit-testable against the closed-form area-mean identity. Its one cost is that it does not repeat out-of-range data — hence the multi-cycle tiling in Phase C4.

### D-008 — Neck metrics in Python, not a `surfaceFieldValue` function object
*2026-08-04*
See Phase C3. The clinical neck inflow rate is a *positive-part* integral, `Σ A·max(U·n̂,0)`, which no `surfaceFieldValue` operation can express; `areaNormalIntegrate` yields net flux, which is ~0 through a sealed sac neck by conservation. Fixing the FOs would have produced a correct implementation of the wrong quantity. Slicing the volume field in Python also gives the clipped disc that KEL needs, so the deferred KEL biomarker is now one line away.

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

**Phase C2 is VALIDATED and SHIPPED in v0.1.0 (2026-06-17).** Public release lives on `main`;
ongoing work (and these dev docs) on `development`. See the repo's two-branch topology.

**Phase C3 and C4 are IMPLEMENTED on branch `phase-c` (2026-08-04), not yet merged.**
Both are fully unit-tested but neither has been run against a real solved case.

**Next step: validate on a real case, then merge `phase-c` → `development`.**
1. `--postprocess-only <case> --neck-metrics` on a solved aneurysm case; check
   `validation.net_flux_near_zero` and confirm the disc in ParaView (full procedure in
   the Phase C3 section above).
2. A 3-cycle `--womersley` run: confirm the inlet flux tracks the waveform across
   **all** cycles (this is what the multi-cycle boundaryData fix addresses) and that
   WSS/OSI stay plausible against the parabolic baseline.
3. Flip `# PHASE-C3-VALIDATION-SWITCH` in `cli.py`, then merge.

**Deferred — WSSG + KEL (pvbatch).** See Phase C2 notes and
`docs/planning/vortexcfd_biomarkers_plan.md` §6. C3 leaves KEL one line away: it is
`0.5·ρ·|U|²·(U·n̂)` integrated over the disc `neck.disc_sample` already returns.

**Roadmap order after that: Phase D (extensions) → Phase E (extended biomarker suite,
see `docs/planning/Biomarcadores_candidatos.md`).** Note the geometry-only biomarkers
(SR, NSI, UI) are out of scope here — tracked upstream as `Supperpt/VORTEX` milestone #2.

---

## 8. Session history

| Date | Session summary |
|---|---|
| 2026-05-26 | Phase A implementation: all Python modules + 13 OpenFOAM templates + Allrun script. Syntax-validated. Not yet run against real STLs on Linux. |
| 2026-05-26 | pytest suite: 99 tests across 4 files (waveform, scaling, case_builder, env_check). All pass on Windows with synthetic pyvista geometry. Added pyproject.toml. |
| 2026-05-29 | Phase A validation (partial): pytest confirmed 99/99 pass on Linux in vortex-aneurysm env. Fixed smoke_test.sh (OUT_DIR/STL_DIR/REPO_ROOT not exported — Python subprocess couldn't read them via os.environ). Smoke test PASSED with real VMTK STLs. OpenFOAM v2406 confirmed at standard path. Full mesher+solver run not yet executed. |
| 2026-05-30 | Phase A fully validated on Kubuntu desktop (Ryzen 5 5600X, OpenFOAM v2406). Fixed 6 bugs during first real run (see Section 5). Key fixes: run-cfd.sh conda/venv detection, `-m vortex_cfd` entry point, background patch in 0/U and 0/p, div(nuEff) in fvSchemes, snappyHexMesh serial-only workaround for v2406 segfault, locationInMesh replaced with inlet-centroid method. Velocity field confirmed inside vessel lumen in ParaView. **Phase A COMPLETE.** |
| 2026-06-04 | Phase C implemented and fully validated on a real patient geometry (OpenFOAM v2406, 6 cores, 3 cycles). `metrics_report.json`: TAWSS mean 9.83 Pa, OSI mean 0.019, low-WSS area 0.17%, high-OSI area 0.84% — all clinically plausible. `wallShearStress` confirmed visible in ParaView with correct pulsatile temporal behaviour. Also fixed `run-cfd.sh` conda env detection (base env / space-in-path bugs). **Phase C COMPLETE.** |
| 2026-06-06 | Phase C2 implemented (branch `aneurysm_dome_biomarkers`): two-patch wall mode (aneurysm_sac + parent_vessel), `--legacy-no-aneurysm` fallback flag, neck-plane function objects, normalised WSS, sac pressure metrics, neck flow rate. 8 source files + 3 test files updated. New unit tests (pure numpy/pyvista, no OpenFOAM). Pending real-case validation with updated VORTEX output. WSSG + KEL deferred to pvbatch script (see `docs/planning/vortexcfd_biomarkers_plan.md` §6). |
| 2026-06-10 | Phase C2 first real run on updated VORTEX output. Fixed BUG-007 (`surfaceFieldValue` FOs used invalid ESI v2406 syntax — `surfaceType`/`patches`/missing `writeFields`/bad plane form — verified correct keys against installed v2406 caseDicts + source; neck flux `sum`→`areaNormalIntegrate`) and BUG-008 (parser read `surface_fieldValue.dat`, actual is `surfaceFieldValue.dat`). Also: `case_builder` now accepts `output_neck_plane.json` and scales the neck origin mm→m; added `from __future__ import annotations` across modules for Python 3.9 (vortex-aneurysm conda env); `requires-python >=3.9`. Pipeline ran to completion through the solver. Test suite 139/140 (the 1 failure is the pre-existing pyvista `locationInMesh` fixture artefact). Remaining: eyeball `metrics_report.json` values + ParaView patches. |
| 2026-06-13 | Restored BUG-007/008 + Python-3.9 LLM.md docs that the Ford commit (`ab6a2f6`, made on a stale checkout) had reverted — the *code* fixes were never lost, only the documentation. Added default STL naming scheme + auto-labelling for unattended batch runs (`patch_labeller._auto_label`): `aneurysm.stl`/`wall.stl`/`inlet.stl`/`outlet_<N>.stl` auto-label with no prompt when all files match; any mismatch falls back to prompting all files (D-006). 11 new tests in `tests/test_patch_labeller.py`; suite 151 pass + 1 pre-existing fixture failure. |
| 2026-06-13 (pm) | Batch of 3 patient runs (AA_010, AA_004, AA_011) at 0.30 m/s ICA inlet, 3 cycles, 6 cores, postprocess. AA_010/AA_004 completed; AA_011 false-aborted at checkMesh (BUG-009: skewness 6.66 > old gate of 4, though mesh was fine). Researched mesh-quality thresholds (OpenFOAM guide + CFD Support: skewness ≤ 20 usable; aneurysm-CFD lit PMC5469453). Raised skewness gate 4 → 20 (= OpenFOAM `maxBoundarySkewness`) and fixed the latent non-ortho regex that never matched checkMesh output. Added `tests/test_runner.py` (5 tests). Suite 156 pass + 1 pre-existing fixture failure. AA_011 not yet re-run (user will do it). |
| 2026-06-14 | Re-ran AA_011 (MRI-derived) → pimpleFoam GAMG pressure FPE on first iteration (BUG-013; run_log.md captured it). Diagnosed as rough-surface → distressed BL (74% sac coverage) → sliver cells → GAMG singular. Decided NOT to fix the solver now (sampling phase); root cause handled upstream in VORTEX remesh. Filed GitHub tracking: VORTEX milestone "Conclude remesh functionality" (#3 tune defaults, #4 spurious-openings bug from cap_label, #5 end-to-end validation); vortex-cfd milestone "Phase B — Robustness & inlet physics" + BUG-013 (#5, deferred, PCG/DIC fallback proposed). No code changes. |
| 2026-06-13 (eve) | Sanity-checked AA_010/AA_004 metrics → **Phase C2 COMPLETE** (WSS/OSI/normalised-WSS validated, physiologically plausible). Found 3 neck-metric issues, opened **Phase C3** and logged BUG-010 (`neck_peak_velocity_ms` null — parenthesised-vector parser), BUG-011 (`neck_peak_flow_rate` sign-naive `max`), CAVEAT-012 (infinite neck plane integrates whole vessel cross-section, not sac orifice; ≈ICA throughput) — all DEFERRED, not fixed. Added per-case `run_log.md` diagnostic feature in `runner.py` (incremental markdown of every step + mesh verdict; survives abort) + 2 tests. Suite 158 pass + 1 pre-existing fixture failure. |
| 2026-06-23 | Code review. Fixed BUG-014: `scaling._any_in_mm` mm/m detection used distance-from-origin (`max(abs(bounds))`) instead of bounding-box extent — translation-variant. Reassessed the review's "critical" severity → not a live bug for head-DICOM data (origin-distance and extent give identical verdicts for cerebral-aneurysm scans), fixed as hardening. Confirmed via git that the function was never previously altered. Added regression test (metre mesh offset +1.5 m → not flagged); scaling suite 15 pass. |
| 2026-06-12 | Replaced the synthetic hand-tuned default inlet waveform with the literature-standard **Ford et al. (2005)** ICA archetype. Digitised the paper's Table 2 ICA feature points → periodic cubic spline (pure-numpy generator `data/generate_ica_ford2005.py`) → bundled `data/ica_ford2005.csv` (100 pts, mean=1, peak 1.657 @ t_norm 0.12 matching P1=1.66). `waveform.py` now loads the bundled CSV by default; `--waveform` still overrides. Confirmed `--mean-velocity` is the cycle-averaged velocity (Q_mean = U_mean × A_inlet). Updated cli help text, `pyproject.toml` package-data, and `tests/test_waveform.py` (Ford feature-timing/amplitude assertions). 140 pass, 1 pre-existing fixture failure. |
| 2026-07-02 | Code review (Fable). Confirmed the DeepSeek findings were genuinely addressed (constants module, `cycles>=2` guard, atexit temp cleanup, area-weighted-mean zero guard, non-ortho regex). Fixed two new bugs: **BUG-015** (`--postprocess-only` silently assumed 3 cycles — `--cycles` now defaults to `None`, full run substitutes 3, so standalone post-processing derives the window from the case data) and **BUG-016** (`run-cfd.sh` hardcoded `~/miniconda3`; now discovers the conda base via `conda info --base` + fallback list, fixing this machine's miniforge setup). Added `tests/test_cli.py` (3 tests). Refreshed the stale Python-3.9 note (env is now 3.10, `requires-python >=3.10`). Suite 162 pass + 1 xfail. Follow-up: cleaned up the three minor `postprocess.py` polish items (dead TAWSS lower-bound check, coupled parent-vessel weights, earliest-only surfaceFieldValue dir parsing → merge all restart dirs) + regression test. Suite 163 pass + 1 xfail. |
| 2026-08-04 | **Phase C3 implemented** on branch `phase-c`. Closed BUG-010/011/CAVEAT-012 by redesign rather than repair: deleted the neck `surfaceFieldValue` FOs entirely, because the clinical *neck inflow rate* is a positive-part integral (`Σ A·max(U·n̂,0)`) that no FO operation can express, and `areaNormalIntegrate` measures net flux, which averages to ~0 through a sealed sac neck. New leaf module `vortex_cfd/neck.py` fits the orifice to `aneurysm_sac.stl`'s open boundary loop by SVD (planarity 4.8e-7 measured) and slices the volume `U` snapshots. Key hardening: the normal's **sign is derived geometrically** (oriented toward the sac), never from `neck_plane.json`, whose convention is undocumented and whose normal was never normalised — a flipped normal would have reported outflow as inflow, undetectably. Radius is `r_eff` not `r_max` (biases the disc small; over-inclusion re-creates CAVEAT-012). Disc clipped with `clip_scalar` (−1.3 % area error vs +3.4 % for a cell-centre mask). Added `validation.net_to_inflow_ratio` so CAVEAT-012 is now self-detecting without ParaView. Parser hardened: strips `()` and **warns on malformed rows** instead of `except: pass` — the silence was the worse half of BUG-010. Gated behind `--neck-metrics` (default off, `# PHASE-C3-VALIDATION-SWITCH`) since no solved case was available; default `metrics_report.json` is byte-identical to v1.0.0. Suite 163 → 212 pass + 1 xfail. |
| 2026-08-04 (cont.) | **Phase C4 implemented** on branch `phase-c`. Re-applied the Womersley inlet from the abandoned `phase_b_womersley_inlet` branch onto current code, dropping everything unrelated it bundled (print→logging migration — would have deleted `run_log.md` and the BUG-009 regex fallback — plus `--dry-run`, flow-rate validation, snappy retry). Fixed the prototype's real defect: `timeVaryingMappedFixedValue` has no `outOfBounds repeat`, so its single-cycle boundaryData would have **frozen the inlet** at the cycle-1 end value for every later cycle — silent and plausible-looking. Data is now tiled across `cycles·T` with a terminal sample. `_location_in_mesh` promoted to `_inlet_geometry` (centroid/normal/radius/interior — it already computed all four); `build_case` returns `(case_dir, inlet_params)`. boundaryData written between snappyHexMesh and checkMesh, and **aborts** on failure (an unreadable inlet BC would waste hours of solve). Verified the area-mean identity by radial quadrature: 1.0000 mean, ±2 % across the cycle at 8 harmonics. Cost 301 dirs / 7 MB for 3 cycles × 1200 faces. Added `scipy>=1.11`. Suite 212 → 261 pass + 1 xfail. Neither C3 nor C4 has been run against a real OpenFOAM case — that is the merge gate. |
