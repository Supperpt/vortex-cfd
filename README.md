# vortex-cfd

**Automated computational fluid dynamics for cerebral aneurysms — from a patient STL to hemodynamic biomarkers, in one command.**

---

## What this is, and why it exists

Cerebral aneurysms are localised dilatations of intracranial arteries. Whether one will rupture is poorly predicted by size alone — local **hemodynamics** (the way blood flows inside the sac) play a central role. The clinically relevant biomarkers are:

- **WSS** (Wall Shear Stress) — the tangential force the blood exerts on the vessel wall. Both abnormally low *and* abnormally high WSS are associated with rupture, by different mechanisms.
- **OSI** (Oscillatory Shear Index) — measures how much the WSS direction reverses during a cardiac cycle. High OSI marks zones of disturbed, recirculating flow.
- **TAWSS** (Time-Averaged WSS) — the mean WSS over a full cycle, used together with OSI to map "low-and-oscillatory" risk regions.

Computing these biomarkers for a patient requires a **patient-specific CFD simulation**: take the 3D geometry from an Angio-CT, mesh it, solve pulsatile Navier–Stokes, and post-process. The end-to-end workflow is well-established academically, but in practice it involves dozens of fragile manual steps — meshing parameters, boundary conditions, solver tuning, post-processing scripts. Most labs maintain a bespoke pipeline that takes hours of operator time per patient.

**vortex-cfd is that pipeline, automated.** The user provides only what is patient-specific (the geometry, the mean inlet velocity, optionally a measured waveform); everything else — mesh strategy, boundary conditions, solver settings, post-processing — is fixed by sensible, reproducible defaults agreed with the user (see *Numerical choices* below).

This is the second stage of a PhD-research workflow:

```
DICOM (Angio-CT)  ──[VORTEX]──►  Watertight STL  ──[vortex-cfd]──►  WSS / OSI / TAWSS report
```

[VORTEX](../VORTEX) handles segmentation, surface extraction, and capping. This project takes its output and runs the actual fluid simulation.

---

## What the program does, step by step

Given a directory of STLs produced by VORTEX with `--split-patches` (one file per surface — the lumen wall and each capped opening), the program:

1. **Detects and sources the OpenFOAM environment** (target: ESI **v2512** from openfoam.com; accepts v2406/v2412/v2506/v2512).
2. **Interactively asks the user to label each STL** as `wall`, `inlet`, or `outlet`. VORTEX's caps are numbered geometrically by VMTK and carry no semantic information — the user is the only entity who knows which opening of the model corresponds to the parent artery (the inlet) and which are distal branches (the outlets). The labelling is a 30-second step done once per patient.
3. **Scales the geometry from millimetres to metres.** Medical imaging works in mm; OpenFOAM assumes SI metres. Failing to scale gives results that look plausible but are off by 10⁹ in velocity — a class of bug already learned the hard way in a prior iteration of this work.
4. **Generates a complete OpenFOAM case directory** from Jinja2 templates: `0/`, `constant/`, `system/` with `controlDict`, `fvSchemes`, `fvSolution`, `snappyHexMeshDict`, `decomposeParDict`, `meshQualityDict`, `surfaceFeatureExtractDict`, `blockMeshDict`, plus the labelled STLs placed in `constant/triSurface/`.
5. **Meshes the lumen** using `snappyHexMesh` with prismatic boundary layers (essential for accurate WSS — first-cell wall-distance must be small enough that the velocity gradient at the wall is resolved, not approximated).
6. **Runs `pimpleFoam`** — a transient, pressure-implicit, incompressible Navier–Stokes solver — for the requested number of cardiac cycles, with the inlet velocity modulated by the cardiac waveform.
7. **(Phase C, planned)** Computes WSS, TAWSS, and OSI on the wall over the last cycle, and writes a `metrics_report.json` with summary statistics.

The output is a standard OpenFOAM case directory that can also be opened directly in ParaView via the `.foam` placeholder file for visual inspection or further analysis.

---

## Inputs

The user provides three things (one optional):

| Input | Purpose |
|---|---|
| `--stl-dir <dir>` | Directory containing the wall + cap STLs from VORTEX |
| `--cycles <N>` | Number of cardiac cycles to simulate (typical: 3 — the first is discarded as transient, the last is analysed) |
| `--mean-velocity <m/s>` | Average blood velocity at the inlet over a cardiac cycle (typical for ICA: 0.3–0.5 m/s; can come from a 4D-flow MRI measurement or literature) |
| `--waveform <csv>` *(opt.)* | Two-column CSV (`time_normalised`, `flow_normalised`) defining the pulse shape. Default: a Womersley-like analytical shape with systolic peak at ~30% of the cycle |
| `--cores <N>` | CPU cores for parallel meshing and solving (default: all available) |

Everything else — fluid density, viscosity, turbulence model, time-stepping strategy, mesh refinement levels, solver tolerances — is fixed by the *Numerical choices* below.

---

## Numerical choices and why

These are the deliberate, locked decisions for this pipeline. They are not user-tunable in normal operation; deviating from them requires editing the templates and is documented as a research-mode change.

### Geometry and units
- **Scaling:** auto-detected ×0.001 (mm→m) when the bounding-box maximum exceeds 1.0.
- **Inlet/outlet identification:** interactive prompt at runtime. VMTK assigns cap IDs geometrically; only a human knows which opening is the parent vessel.

### Fluid model
- **Density:** ρ = 1060 kg/m³ (whole blood).
- **Viscosity (default):** Newtonian, μ = 0.0035 Pa·s (ν = 3.3 × 10⁻⁶ m²/s). Justified by the relatively high shear rates in the parent artery; the simplification is the academic norm and is the baseline against which non-Newtonian corrections are measured.
- **Viscosity (opt-in for Phase D):** Carreau non-Newtonian model with Cho & Kensey (1991) parameters: μ₀ = 0.056, μ∞ = 0.0035 Pa·s, λ = 3.313 s, n = 0.3568. Relevant in the aneurysm sac where shear rates drop and viscosity rises.

### Flow regime and turbulence
- **Laminar** by default. Reynolds number in the internal carotid artery is typically 200–400 — well below the transition threshold for a straight tube. Disturbed flow in the aneurysm sac may transition locally, but pulsatile-laminar is the established starting point.
- **k-ω SST** available via opt-in flag for cases where transitional/turbulent behaviour is suspected.

### Solver: `pimpleFoam`
- Transient, incompressible Navier–Stokes with the PIMPLE algorithm (combines PISO and SIMPLE — robust for moderate-Δt simulations with sub-iteration coupling between pressure and velocity).
- Rigid-wall (no fluid–structure interaction in this project for now).
- **Time-stepping:** adaptive, driven by `maxCo = 0.8` (the Courant number — a dimensionless measure of how far fluid moves per timestep relative to mesh cell size). Adaptive timestepping is mandatory because peak systolic velocity is ~10× the diastolic baseline; a fixed Δt either wastes effort during diastole or diverges during systole.
- **Schemes:** backward second-order in time, Gauss linearUpwind for advection (stable and ~second-order accurate).

### Boundary conditions
- **Inlet velocity:** `flowRateInletVelocity` for Phase A — applies a uniform parabolic (Poiseuille) profile scaled by the cardiac waveform. Adequate because VORTEX's flow extensions (typically 5× the local radius) let the profile re-develop before reaching the aneurysm.
- **Inlet velocity (Phase B):** Womersley profile available via opt-in flag — implemented via `codedFixedValue` with Fourier decomposition of the waveform and Bessel functions for the radial profile. Required for academic publication where Womersley is the standard.
- **Inlet pressure:** zero-gradient.
- **Outlet velocity:** `inletOutlet` (acts as zero-gradient when flow is outgoing, prevents inflow if recirculation reaches the outlet — a common numerical instability in vascular CFD).
- **Outlet pressure:** fixed at 0 Pa (gauge — only relative pressure matters for incompressible flow).
- **Wall:** no-slip on velocity, zero-gradient on pressure.

### Cardiac cycle
- **Period:** T = 0.857 s (70 bpm) by default.
- **Cycles simulated:** user-specified N; the first is discarded as transient initialisation, the last is used for biomarker averaging.

### Meshing strategy
- `blockMesh` background grid sized automatically from the wall bounding box + 20% buffer.
- `snappyHexMesh` castellated → snap → addLayers, with 4 prismatic boundary layers, expansion ratio 1.3, final-layer-thickness 0.3 of the bulk cell size. Layers are non-negotiable for WSS accuracy — without them, the wall gradient is interpolated across a single large cell and WSS is meaningless.
- `checkMesh` aborts the run if `maxNonOrtho > 70°` or `maxSkewness > 4`. These thresholds catch the kind of pathological cells that crash `pimpleFoam` later.

### Parallelisation
- `decomposePar` with the scotch method.
- All meshing and solving steps run under `mpirun -np <cores>`.

---

## Roadmap

**Phase A — MVP (in progress).** End-to-end from STLs to a runnable, openable OpenFOAM case. Skeleton, CLI, interactive labelling, scaling, Jinja2 templates, mesh, solve, no post-processing. *Success criterion:* a real patient STL goes in, a `case_XXX.foam` comes out that opens in ParaView and shows reasonable velocity fields.

**Phase B — Robustness and Womersley.** Womersley inlet profile (publication-grade), retry logic if `snappyHexMesh` fails, structured logging, validation that the inlet area / mean velocity combination is physiologically plausible.

**Phase C — Post-processing and reports.** WSS, OSI, TAWSS as OpenFOAM function objects with `fieldAverage` over the last cycle. JSON report with summary statistics. Optional ParaView screenshots via `pvbatch`.

**Phase D — Extensions.** Carreau non-Newtonian opt-in. Mesh-independence helper (runs three densities and reports convergence). Folder hooks reserved for future FSI work.

---

## Status

Phase A in progress — repository scaffolding only. No runnable code yet.

---

## Requirements

- Linux (native; not WSL).
- OpenFOAM **v2512** from openfoam.com. *This is the ESI/OpenCFD release, not the Foundation v12 from openfoam.org.* The two are parallel forks of the same project with slightly different dict syntax and function object catalogues. All templates here are written for ESI v2512; v2406, v2412, and v2506 are also accepted.
- Python 3.10+.

## Setup (one-time)

Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html), then:

```bash
conda create -n vortex-aneurysm python=3.12 -y
conda activate vortex-aneurysm
cd vortex-cfd
bash setup.sh
```

## Usage

Every new terminal session, activate the environment first:

```bash
conda activate vortex-aneurysm
```

Then run the pipeline:

```bash
bash run-cfd.sh --stl-dir <path/to/stls> --cycles 3 --mean-velocity 0.4 --cores 4 --out-dir <output-dir>
```

---

## Project layout (planned)

```
vortex-cfd/
├── README.md
├── requirements.txt
├── setup.sh                  Adds Python deps to the vortex-aneurysm conda env
├── run-cfd.sh                Launcher — sources OpenFOAM bashrc, then runs the CLI
├── vortex_cfd/
│   ├── cli.py                Click-based entry point
│   ├── env_check.py          Detects and sources OpenFOAM, validates version
│   ├── patch_labeller.py     Interactive wall/inlet/outlet labelling
│   ├── scaling.py            mm→m detection and STL rewriting
│   ├── waveform.py           Default and user-CSV cardiac waveforms
│   ├── case_builder.py       Renders Jinja2 templates into a case directory
│   ├── runner.py             Orchestrates surfaceFeatureExtract → blockMesh → snappyHexMesh → pimpleFoam → reconstructPar
│   └── templates/            Jinja2 templates for all OpenFOAM dicts
└── tests/
    └── smoke_test.sh         End-to-end run on a known STL
```

---

## Relationship to the prior iteration

A previous attempt at this same goal existed as `VesselForge_AutoCFD`. It is not preserved as code in this repository, but its lessons are. The relevant ones, documented in `../VORTEX/previous_cfd_iteration.md`, are:

- mm→m scaling must be applied to STL geometry *and* to any patch-coordinate metadata files. Forgetting the latter silently maps the inlet to a point 1000× away from the mesh.
- `patch_labels.json` must be written after the case directory exists, not before. The prior race condition produced runs that "completed" against an empty directory.
- Hardcoding 4 CPU cores leaves modern hardware idle. Expose `--cores` from the start.
- Validation ranges from the prior work, used here as sanity checks: WSS 0–50 Pa, OSI 0–0.5, velocity 0.1–1.0 m/s, pressure 0–200 Pa relative to the outlet.
