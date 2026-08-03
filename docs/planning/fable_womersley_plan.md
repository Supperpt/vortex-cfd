# Phase B — Womersley inlet profile: implementation plan

## Context

Phase B (publication-grade Womersley inlet) was prototyped on the branch
`phase_b_womersley_inlet` but abandoned. That branch forked from `145156f`, which
predates Phase C2 (two-patch aneurysm mode), the 1.0.0 release, `constants.py`,
`run_log.md`, the checkMesh-regex fixes, and BUG-014/015/016. It sits ~18 commits
behind `development`, so a merge/rebase would regress shipped work. This plan
**re-applies the branch's validated Womersley design onto current `development`**,
adapting it to the current code and dropping the unrelated items the branch bundled.

Goal: an opt-in `--womersley` flag that replaces the default parabolic
`flowRateInletVelocity` inlet with the exact Womersley analytical profile, required
for academic publication.

**Confirmed scope:** Womersley only. The other Phase B items the branch bundled
(`--dry-run`, physiological flow validation, snappyHexMesh retry, and the
`print`→`logging` migration) are out of scope and left as separate follow-ups — the
logging migration in particular would regress `development`'s `run_log.md` feature.
**Verification target:** Python unit tests + case-generation smoke assertions; the
real multi-cycle OpenFOAM run is left for the user to execute.

## What the branch got right (keep)

- **`vortex_cfd/womersley.py`** — pure math + OpenFOAM writer, no OpenFOAM/pyvista
  in the math path:
  - `fourier_coefficients(waveform, N_harm=8)` — FFT of the normalised waveform.
  - `_womersley_shape(r_norm, alpha_k)` — complex Bessel shape `1 - J0(Λs)/J0(Λ)`,
    **area-normalised** so each harmonic's area-mean = 1 (`scipy.special.jv` takes
    complex arguments directly; no polynomial approximations).
  - `womersley_velocities(...)` — reconstructs `U(r,t)` = Poiseuille DC term +
    `2·Re[Σ C_k·G_k(r)·e^{ikωt}]`; the area-mean reproduces `U_mean·waveform(t)`.
    Near-wall retrograde flow is deliberately not clipped (physiologically correct).
  - `write_boundary_data(...)` — writes `constant/boundaryData/inlet/{points, <t>/U}`.
- **Delivery via `timeVaryingMappedFixedValue`** reading per-face boundaryData — no
  runtime C++ compilation; the math stays plain, testable Python.
- **Post-mesh timing**: boundaryData is written *after* snappyHexMesh (inlet face
  centres only exist then) and *before* checkMesh.
- **`0/U.j2`**: `{% if womersley %}` switches the inlet block to
  `timeVaryingMappedFixedValue`, else the existing `flowRateInletVelocity`.
- **scipy dependency** (`scipy>=1.11`).

## What to drop / change from the branch

- **Do NOT take the branch's `runner.py` / `cli.py` / `__main__.py` wholesale.**
  They migrate `print`→`logging` and remove the `run_log.md` incremental logging and
  the checkMesh fail-loud regex fallback that `development` added later (BUG-009).
  Re-apply only the Womersley hooks on top of the *current* runner.
- **Drop the out-of-scope bundled items**: `--dry-run`, physiological flow-rate
  validation, snappyHexMesh retry (`_run_snappy`), and the logging migration.
- **Two-patch mode**: the branch's `0/U.j2` used a single `{{ wall_patch }}`; the
  current template loops `{% for patch in wall_patches %}`. Keep the loop; only the
  inlet block changes.
- **Reuse current helpers**: read STLs via `scaling.read_stl` (not raw `pv.read`);
  take `NU`/`RHO` from `vortex_cfd/constants.py` (not a hardcoded `3.3e-6`); locate
  the meshed inlet patch with `postprocess._wall_block(mb, "inlet")` instead of an
  ad-hoc multiblock loop.
- **Keep `_inlet_area`** (actual STL area + finite/positive guard) for the parabolic
  flow table; use the equivalent radius `sqrt(A/π)` only for the Womersley profile.

## Correctness fix to add (branch bug)

`timeVaryingMappedFixedValue` does **not** repeat out-of-range data (there is no
equivalent of `flowRateInletVelocity`'s `outOfBounds repeat`). The branch wrote
boundaryData for a single cycle `[0, T]`, so `--cycles > 1` would freeze the inlet
at the cycle-1 end value for every later cycle. Fix: generate/tile the boundaryData
across the **full** `end_time = cycles · T` (replicate the normalised cycle `cycles`
times with shifted absolute times), or make `womersley_velocities` accept
`cycles`/`end_time` and emit the repeated series. Add a regression test asserting the
last written boundaryData time `>= end_time`.

## Files to change

- `vortex_cfd/womersley.py` — **new**, adapted from the branch (constants import,
  multi-cycle tiling).
- `vortex_cfd/templates/0/U.j2` — add the `{% if womersley %}` inlet branch; keep the
  `wall_patches` loop.
- `vortex_cfd/case_builder.py` — refactor `_location_in_mesh` → `_inlet_geometry`
  (returns centroid, inward normal, radius, interior point) with `_location_in_mesh`
  kept as a thin wrapper for existing callers/tests; add a `womersley` param; add
  `"womersley"` to the Jinja context; return `(case_dir, inlet_params)` where
  `inlet_params` carries centroid / normal / radius / mean_velocity / waveform / nu
  for the runner.
- `vortex_cfd/runner.py` — `run_pipeline(..., womersley=False, inlet_params=None)`;
  after the existing serial `snappyHexMesh` step and before `checkMesh`, when
  `womersley`, read inlet face centres from the mesh and call
  `womersley.womersley_velocities` + `write_boundary_data`. Leave `run_log.md`
  logging and the mesh-quality gate intact.
- `vortex_cfd/cli.py` — add the `--womersley` flag; unpack the new `build_case`
  tuple; forward `womersley` / `inlet_params` to `run_pipeline`. (No logging change.)
- `requirements.txt` / `pyproject.toml` — add `scipy>=1.11`.
- `tests/test_womersley.py` — port the branch's tests (Fourier coefficients, shape
  function, velocities, boundaryData writer) + a new multi-cycle-coverage test.
- `tests/test_case_builder.py` — update for the `build_case` tuple return and
  `_inlet_geometry`.
- `LLM.md` — Phase B section: mark Womersley implemented; record the design and the
  multi-cycle fix.

## Alternative considered — `codedFixedValue` (rejected)

The original README / LLM.md wording implied a `codedFixedValue` inlet (Fourier +
Bessel evaluated inside an OpenFOAM coded BC). It avoids writing boundaryData files
and auto-covers all simulation time. Rejected because: it requires runtime C++
compilation (`dynamicCode`, fragile across ESI versions and on locked-down HPC
nodes); OpenFOAM has no native complex Bessel, so `J0` of a complex argument would
have to be hand-coded in C++; and the math becomes untestable in Python. The
branch's `timeVaryingMappedFixedValue` approach is compilation-free, unit-testable,
and already validated — the better fit, and almost certainly what made the original
attempt feel "tricky."

## Verification

1. `pytest` — full suite green, including new `test_womersley.py` (math + writer +
   multi-cycle coverage) and the updated `test_case_builder.py`.
2. Case-generation smoke: build a case with `--womersley` on a sample STL dir and
   assert `0/U` renders `timeVaryingMappedFixedValue`; after a mesh, assert
   `constant/boundaryData/inlet/` contains `points` + per-time `U` directories
   spanning `end_time`. A `--womersley` vs. non-womersley diff of `0/U` confirms the
   switch.
3. (User-run, optional) Real OpenFOAM run: a 3-cycle `--womersley` case reaches the
   solver, the inlet flux tracks the waveform across *all* cycles (confirming the
   multi-cycle fix), and WSS/OSI stay physiologically plausible versus the parabolic
   baseline.
