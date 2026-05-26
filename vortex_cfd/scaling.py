"""Detect millimetre geometry and scale all STLs to metres (SI)."""

import tempfile
from pathlib import Path

import pyvista as pv

MM_TO_M = 0.001
# If any axis of the bounding box exceeds this value (metres assumed), it's in mm.
_BBOX_THRESHOLD = 1.0


def _any_in_mm(stl_paths: list[Path]) -> bool:
    for path in stl_paths:
        mesh = pv.read(str(path))
        if max(abs(v) for v in mesh.bounds) > _BBOX_THRESHOLD:
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
