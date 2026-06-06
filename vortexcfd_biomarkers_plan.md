# vortex-cfd Plan: Aneurysm-Scoped Biomarker Computation

## Goal

Extend the pipeline to accept `aneurysm_sac.stl` + `parent_vessel.stl` (+ `neck_plane.json`) from VORTEX, mesh them as named patches, and compute all PhD biomarkers scoped exclusively to the aneurysm sac surface. No new solver physics — this is a post-processing and meshing labelling change.

## Prerequisite

VORTEX must be updated first (see `VMTK_plan_biomarkers.md`) to produce:
- `aneurysm_sac.stl` — aneurysm dome, clipped at the neck
- `parent_vessel.stl` — remainder of the lumen wall (replaces `wall.stl`)
- `neck_plane.json` — `{"origin": [...], "normal": [...]}` at the neck cross-section

---

## 1. patch_labeller.py — Accept Two Wall Patches

Currently the labeller enforces exactly **1 wall + 1 inlet + ≥1 outlet**. The new contract is:
- **1 `aneurysm_sac`** (wall type)
- **1 `parent_vessel`** (wall type)  
- **1 inlet**
- **≥1 outlet**

`patch_labeller.py` should present `aneurysm_sac` and `parent_vessel` as distinct options with the `wall` role, and write both into `patch_labels.json` with the correct keys.

---

## 2. snappyHexMeshDict.j2 — Register aneurysm_sac as a Named Patch

In the `geometry` block (already iterates `all_stls`): no change needed if `aneurysm_sac` is included in `scaled_stls`.

In `refinementSurfaces`, add an explicit entry for `aneurysm_sac` with the same wall treatment as the existing wall patch:

```
aneurysm_sac
{
    level       (3 4);
    patchInfo   { type wall; inGroups (wall); }
}
```

And rename the existing `wall` entry to `parent_vessel` (or keep both names depending on final naming convention from VORTEX).

In `addLayersControls`, add `aneurysm_sac` to the layer-addition surface list alongside `parent_vessel` so both get the 4-layer prismatic BL treatment needed for WSS accuracy.

---

## 3. 0/U.j2 and 0/p.j2 — Add aneurysm_sac BC Entries

`0/U.j2`:
```
aneurysm_sac
{
    type    noSlip;
}
```

`0/p.j2`:
```
aneurysm_sac
{
    type    zeroGradient;
}
```

Same entries for `parent_vessel` (replacing the old `wall` entry).

---

## 4. controlDict.j2 — Extended Function Objects

### 4a. wallShearStress — Compute on Both Wall Patches

```
wallShearStress
{
    type            wallShearStress;
    libs            (fieldFunctionObjects);
    patches         (aneurysm_sac parent_vessel);
    writeControl    writeTime;
}
```

### 4b. surfaceFieldValue — Neck Inflow Rate and Peak Velocity

The neck plane geometry comes from `neck_plane.json` (passed into the Jinja2 context as `neck_origin` and `neck_normal`).

```
surfaceFieldValue_neck_flux
{
    type            surfaceFieldValue;
    libs            (fieldFunctionObjects);
    surfaceType     plane;
    basePoint       ( {{ neck_origin[0] }} {{ neck_origin[1] }} {{ neck_origin[2] }} );
    normalVector    ( {{ neck_normal[0] }} {{ neck_normal[1] }} {{ neck_normal[2] }} );
    fields          (U);
    operation       sum;
    writeControl    writeTime;
}

surfaceFieldValue_neck_peak_vel
{
    type            surfaceFieldValue;
    libs            (fieldFunctionObjects);
    surfaceType     plane;
    basePoint       ( {{ neck_origin[0] }} {{ neck_origin[1] }} {{ neck_origin[2] }} );
    normalVector    ( {{ neck_normal[0] }} {{ neck_normal[1] }} {{ neck_normal[2] }} );
    fields          (U);
    operation       max;
    writeControl    writeTime;
}
```

These write time-series CSVs to `postProcessing/` automatically.

### 4c. surfaceFieldValue — Sac Pressure (Mean and Peak)

```
surfaceFieldValue_sac_pressure_mean
{
    type            surfaceFieldValue;
    libs            (fieldFunctionObjects);
    surfaceType     patch;
    patches         (aneurysm_sac);
    fields          (p);
    operation       areaAverage;
    writeControl    writeTime;
}

surfaceFieldValue_sac_pressure_max
{
    type            surfaceFieldValue;
    libs            (fieldFunctionObjects);
    surfaceType     patch;
    patches         (aneurysm_sac);
    fields          (p);
    operation       max;
    writeControl    writeTime;
}
```

Physical pressure Pa = `p_kinematic × rho` (same conversion as WSS, already in `postprocess.py`).

---

## 5. postprocess.py — Scope Metrics to aneurysm_sac Patch

`read_wss_series()` currently extracts the `wall` block from the OpenFOAM reader. Update to:
- Read `aneurysm_sac` block for sac-scoped TAWSS/OSI.
- Read `parent_vessel` block for parent-vessel mean WSS (used as denominator for normalised WSS).
- Normalised WSS = `TAWSS_sac_local / TAWSS_parent_mean`.

`compute_metrics()` should also:
- Parse the `postProcessing/surfaceFieldValue_neck_flux/` time-series → neck inflow rate Q(t), peak Q.
- Parse `postProcessing/surfaceFieldValue_sac_pressure_*/` → mean/peak pressure.
- Write all new values into `metrics_report.json` under new keys.

`case_builder.py` needs to:
- Read `neck_plane.json` from the STL directory and pass `neck_origin` + `neck_normal` into the Jinja2 template context.

---

## 6. WSSG and KEL — ParaView pvbatch (Deferred)

These two metrics require surface spatial gradients and volumetric flux integration; there is no native OpenFOAM function object for either.

**Standard method (pvbatch script):**

WSSG:
1. Open `.foam` → Extract Block → `aneurysm_sac`.
2. Apply "Gradient Of Unstructured Dataset" on `wallShearStress`.
3. Calculator: `mag(gradient)` → WSSG field.
4. Temporal Statistics → time-averaged WSSG.

KEL:
1. Slice at the neck plane (use `neck_plane.json` coordinates).
2. Calculator: `0.5 * 1060 * mag(U)^2 * dot(U, neck_normal_vec)`.
3. Integrate Variables → KEL(t) time series.
4. Report cycle-averaged value.

These are implemented as a `pvbatch` Python script called from `runner.py` after the solve, gated by a new `--screenshots` / `--pvbatch` flag (already deferred from Phase C).

---

## 7. case_builder.py — Context Updates

New Jinja2 context variables to add:
```python
"aneurysm_patch": "aneurysm_sac",
"parent_vessel_patch": "parent_vessel",
"neck_origin": [x, y, z],      # from neck_plane.json
"neck_normal": [nx, ny, nz],    # from neck_plane.json
```

Load `neck_plane.json` from the STL source directory in `build_case()`.

---

## Metrics Summary

| Metric | Source | Location |
|---|---|---|
| TAWSS (sac) | postprocess.py via pyvista | `aneurysm_sac` patch |
| Max WSS | postprocess.py | `aneurysm_sac` patch, peak systole |
| Normalised WSS | postprocess.py | sac TAWSS / parent mean TAWSS |
| OSI (sac) | postprocess.py via pyvista | `aneurysm_sac` patch |
| WSSG | pvbatch script | `aneurysm_sac` patch surface gradient |
| Neck inflow rate Q(t) | OpenFOAM `surfaceFieldValue` | neck cutting plane |
| Peak inflow velocity | OpenFOAM `surfaceFieldValue` | neck cutting plane |
| Mean static pressure | OpenFOAM `surfaceFieldValue` | `aneurysm_sac` patch |
| Peak static pressure | OpenFOAM `surfaceFieldValue` | `aneurysm_sac` patch |
| Kinetic Energy Loss | pvbatch script | neck cutting plane integration |

---

## Implementation Order

1. Confirm VMTK side produces the three required files (`aneurysm_sac.stl`, `parent_vessel.stl`, `neck_plane.json`) with the naming convention above.
2. `patch_labeller.py` — accept two wall patches.
3. `snappyHexMeshDict.j2` — add `aneurysm_sac` refinement surface + layers.
4. `0/U.j2`, `0/p.j2` — add `aneurysm_sac` BC entries.
5. `controlDict.j2` — add `surfaceFieldValue` objects (neck plane + sac pressure).
6. `case_builder.py` — load `neck_plane.json`, pass new context variables.
7. `postprocess.py` — scope TAWSS/OSI to sac patch, parse new `postProcessing/` CSVs, update `metrics_report.json`.
8. pvbatch script for WSSG + KEL (separate, deferred).
