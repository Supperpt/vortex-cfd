"""Detect millimetre geometry and scale all STLs to metres (SI)."""
from __future__ import annotations

import tempfile
from pathlib import Path

import pyvista as pv

MM_TO_M = 0.001
# If any bounding-box extent exceeds this value (metres assumed), it's in mm.
_BBOX_THRESHOLD = 1.0


def _any_in_mm(stl_paths: list[Path]) -> bool:
    for path in stl_paths:
        mesh = pv.read(str(path))
        # Use bounding-box extents (size), not max absolute coordinate, so the
        # check is translation-invariant. Medical scans carry coordinates tied to
        # the scanner isocenter, so measuring distance-from-origin could misjudge
        # the unit if the geometry is offset from (0,0,0).
        xmin, xmax, ymin, ymax, zmin, zmax = mesh.bounds
        max_extent = max(xmax - xmin, ymax - ymin, zmax - zmin)
        if max_extent > _BBOX_THRESHOLD:
            return True
    return False


def scale_stls(stl_paths: list[Path], labels: dict[Path, str]) -> dict[str, Path]:
    """
    Copy all STLs to a temp directory under canonical names (wall / inlet / outlet_N).
    Scale ×0.001 when millimetre coordinates are detected.

    Returns {canonical_name: Path} — e.g. {'wall': ..., 'inlet': ..., 'outlet_0': ...}.

    The mm→m scaling must be applied to geometry *and* written back to file so that
    snappyHexMesh and the boundary-condition patches are consistent in SI metres.
    """
    needs_scale = _any_in_mm(stl_paths)
    if needs_scale:
        print("Detected millimetre coordinates — scaling ×0.001 (mm → m).")
    else:
        print("Coordinates appear to be in metres — no scaling applied.")

    tmp = Path(tempfile.mkdtemp(prefix="vortex_cfd_stls_"))
    outlet_idx = 0
    result: dict[str, Path] = {}

    for path in stl_paths:
        label = labels[path]
        if label == "outlet":
            canonical = f"outlet_{outlet_idx}"
            outlet_idx += 1
        else:
            canonical = label  # 'wall' or 'inlet'

        mesh = pv.read(str(path))
        if needs_scale:
            mesh.points *= MM_TO_M

        dest = tmp / f"{canonical}.stl"
        mesh.save(str(dest), binary=False)
        result[canonical] = dest

    return result
