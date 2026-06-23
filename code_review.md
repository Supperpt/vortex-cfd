# `vortex-cfd` Code Review Report

## 1. Overview and Architecture
`vortex-cfd` is an elegant, modular pipeline for automating OpenFOAM CFD simulations for cerebral aneurysms. The project leverages `pyvista` for mesh manipulation, `Jinja2` for generating OpenFOAM dictionaries, and `click` for a clean CLI.

The architecture is well-separated into logical phases:
- **Environment & Pre-processing:** `env_check.py`, `scaling.py`, `patch_labeller.py`
- **Case Generation:** `case_builder.py`, `waveform.py`, Jinja2 templates
- **Execution Orchestration:** `runner.py`
- **Post-processing & Biomarkers:** `postprocess.py`

---

## 2. Strengths and Good Practices

### Robust Separation of Concerns
The project perfectly isolates IO operations from mathematical operations. In `postprocess.py`, the pure math functions (`tawss`, `osi`, `summary_stats`) depend solely on Numpy arrays and not on PyVista or OpenFOAM readers. This makes unit testing trivial and ensures long-term maintainability.

### Clever Automation & Heuristics
- **Inlet Center Detection:** In `case_builder.py`, `_location_in_mesh` correctly calculates a reliable point strictly inside the fluid domain by stepping inwards from the inlet's area-weighted centroid along its negative normal. This brilliantly avoids the pitfalls of bounding-box centers falling outside the fluid domain in curved vessels.
- **Waveform Normalization:** `waveform.py` forces the integrated mean of the waveform to `1.0`. This guarantees that the user-provided `--mean-velocity` dictates the exact physical flow rate, irrespective of the underlying shape's native scale.

### OpenFOAM Integration & Safety
- **Quality Gate:** `runner.py` uses `checkMesh` output to parse `maxNonOrtho` and `maxSkewness`. Aborting execution when limits are exceeded (>70° and >20) saves time by preventing inevitable `pimpleFoam` crashes.
- **Defensive Meshing:** The comments indicate that `snappyHexMesh` is intentionally run in serial to avoid a known load-balancing segmentation fault bug in OpenFOAM ESI v2406. This demonstrates excellent domain knowledge and defensive programming.
- **Environment Detection:** `env_check.py` seamlessly searches for OpenFOAM in multiple directories (`/opt`, `/usr/lib/openfoam`) and accurately captures variables from `bashrc` using a subprocess sub-shell.

---

## 3. Critical Issues & Bugs

### 🚨 **CRITICAL: Coordinate Translation Bug in `scaling.py`**
The heuristic used to detect if STLs are in millimeters vs. meters is flawed and will produce catastrophic scaling errors for certain patient datasets.
In `scaling.py`:
```python
def _any_in_mm(stl_paths: list[Path]) -> bool:
    for path in stl_paths:
        mesh = pv.read(str(path))
        if max(abs(v) for v in mesh.bounds) > _BBOX_THRESHOLD: # Threshold is 1.0
            return True
    return False
```
**The Problem:** `mesh.bounds` returns `(xmin, xmax, ymin, ymax, zmin, zmax)`. By checking `max(abs(v))`, you are measuring the **absolute distance from the origin**, not the size of the mesh. In medical imaging (MRI/CT), coordinates are tied to the scanner's isocenter. A brain mesh might be correctly sized in meters (e.g., 0.2m wide), but positioned at `Z = 1.5m` from the scanner origin. Your logic would falsely flag this mesh as being in "millimeters" (since 1.5 > 1.0) and erroneously scale it down by 1000.

**The Fix:** You should evaluate the **extents** (dimensions) of the bounding box, making the check translation-invariant:
```python
        bounds = mesh.bounds
        max_extent = max(bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4])
        if max_extent > _BBOX_THRESHOLD:
```

### ⚠️ **Warnings on Error Handling in Subprocesses**
In `runner.py`, `subprocess.run(cmd, shell=True, ...)` is heavily used.
While `cores` and `cycles` are validated by `click` as integers (preventing standard shell injections), using `shell=True` with dynamically assembled strings is generally considered a minor security risk and prone to quote-escaping bugs if directory paths (like `case_dir`) contain spaces. Since `case_dir` is generated with a timestamp it shouldn't contain spaces natively, but `out_dir` provided by the user might.
**Recommendation:** Enclose paths in quotes when passed into shell strings, or refactor to avoid `shell=True` if possible.

---

## 4. Minor Improvements and Polish

1. **Jinja2 Template Permissions:** In `case_builder.py`, `Allrun.j2` is copied and later given executable permissions using `.chmod()`. This works, but relying on Python to manipulate file permissions across different environments can sometimes have edge cases. It is fine as written, just something to keep an eye on.
2. **Commented-Out Neck Metrics:** In `postprocess.py`, lines 432-443 disable the neck metrics. The documentation indicates this is intentional pending validation, but be sure to track this technical debt explicitly via GitHub Issues to prevent it from remaining permanently disabled.
3. **OSI Div-by-Zero Handling:** In `postprocess.py`, you correctly use `np.errstate(divide="ignore", invalid="ignore")` alongside `np.where(mean_mag > 0, ...)`. This is a clean, pythonic way to handle static flow regions. 
4. **Click Path Checking:** In `cli.py`, you use `click.Path(exists=True)` for `stl_dir`. This correctly offloads error handling to the CLI framework, making the terminal output cleaner.

## Conclusion
Overall, the codebase is highly robust, scientifically sound, and demonstrates a deep understanding of OpenFOAM's quirks and Python tooling. The architecture allows for safe restarts and post-processing recalculations. **The only immediate action required is fixing the bounding-box logic in `scaling.py`** to prevent catastrophic scaling errors when handling translated coordinate systems.
