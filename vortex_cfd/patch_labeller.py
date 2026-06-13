"""Wall / inlet / outlet labelling of VORTEX STL caps.

Two paths:
  * auto-labelling from the default naming scheme (enables unattended batch runs);
  * interactive prompting (fallback when the scheme is not fully matched).
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import click
import pyvista as pv

# Legacy mode: single wall patch.
VALID_LABELS_LEGACY = ("wall", "inlet", "outlet")

# New (default) mode: aneurysm sac + parent vessel as separate wall patches.
VALID_LABELS_NEW = ("aneurysm_sac", "parent_vessel", "inlet", "outlet")

# Default naming scheme: input STL filename stem (case-insensitive) -> label.
# An unattended run requires every STL in the directory to match this scheme;
# otherwise the labeller falls back to full interactive prompting.
#   New mode:    aneurysm.stl -> aneurysm_sac, wall.stl -> parent_vessel,
#                inlet.stl -> inlet, outlet_<N>.stl -> outlet
#   Legacy mode: wall.stl -> wall, inlet.stl -> inlet, outlet_<N>.stl -> outlet
_NAMING_SCHEME_NEW = {"aneurysm": "aneurysm_sac", "wall": "parent_vessel", "inlet": "inlet"}
_NAMING_SCHEME_LEGACY = {"wall": "wall", "inlet": "inlet"}
# Outlets: outlet, outlet_1, outlet_2, ... (the trailing index is optional).
_OUTLET_RE = re.compile(r"^outlet(_\d+)?$", re.IGNORECASE)


def _label_from_filename(path: Path, legacy: bool) -> str | None:
    """Map a single STL filename to a label via the default scheme, or None."""
    stem = path.stem.lower()
    if _OUTLET_RE.match(stem):
        return "outlet"
    scheme = _NAMING_SCHEME_LEGACY if legacy else _NAMING_SCHEME_NEW
    return scheme.get(stem)


def _validate_labels(labels: dict[Path, str], legacy: bool) -> list[str]:
    """Return a list of validation error strings (empty == valid)."""
    errors: list[str] = []
    if legacy:
        walls = [p for p, l in labels.items() if l == "wall"]
        if len(walls) != 1:
            errors.append(f"Expected exactly 1 wall, got {len(walls)}: {[p.name for p in walls]}")
    else:
        sacs = [p for p, l in labels.items() if l == "aneurysm_sac"]
        vessels = [p for p, l in labels.items() if l == "parent_vessel"]
        if len(sacs) != 1:
            errors.append(f"Expected exactly 1 aneurysm_sac, got {len(sacs)}: {[p.name for p in sacs]}")
        if len(vessels) != 1:
            errors.append(f"Expected exactly 1 parent_vessel, got {len(vessels)}: {[p.name for p in vessels]}")

    inlets = [p for p, l in labels.items() if l == "inlet"]
    outlets = [p for p, l in labels.items() if l == "outlet"]
    if len(inlets) != 1:
        errors.append(f"Expected exactly 1 inlet, got {len(inlets)}: {[p.name for p in inlets]}")
    if not outlets:
        errors.append("Expected at least 1 outlet.")
    return errors


def _auto_label(stl_paths: list[Path], legacy: bool) -> dict[Path, str] | None:
    """
    Label every STL from the default naming scheme.

    Returns the {Path: label} dict only if *every* file matches the scheme and the
    result passes validation; otherwise returns None so the caller prompts instead.
    """
    labels: dict[Path, str] = {}
    for path in stl_paths:
        label = _label_from_filename(path, legacy)
        if label is None:
            return None  # unrecognised filename -> not a clean scheme match
        labels[path] = label
    if _validate_labels(labels, legacy):
        return None  # matched names but wrong counts -> fall back to prompting
    return labels


def _describe(path: Path) -> str:
    try:
        mesh = pv.read(str(path))
        c = mesh.center
        sized = mesh.compute_cell_sizes()
        area = float(sized.cell_data["Area"].sum())
        b = mesh.bounds
        return (
            f"  centre ({c[0]:.2f}, {c[1]:.2f}, {c[2]:.2f})  "
            f"area {area:.2f}  "
            f"bbox x[{b[0]:.2f},{b[1]:.2f}] y[{b[2]:.2f},{b[3]:.2f}] z[{b[4]:.2f},{b[5]:.2f}]"
        )
    except Exception:
        return "  (could not read geometry)"


def label_patches(stl_paths: list[Path], legacy: bool = False) -> dict[Path, str]:
    """
    Prompt the user to assign a label to each STL.

    New mode (default): labels are aneurysm_sac / parent_vessel / inlet / outlet.
      Validates: exactly 1 aneurysm_sac + 1 parent_vessel + 1 inlet + ≥1 outlet.

    Legacy mode (--legacy-no-aneurysm): labels are wall / inlet / outlet.
      Validates: exactly 1 wall + 1 inlet + ≥1 outlet.

    Returns a {Path: label} dict.

    If every STL matches the default naming scheme (see module docstring), labels
    are assigned automatically with no prompting — enabling unattended batch runs.
    Any unrecognised filename or count mismatch falls back to prompting for all files.
    """
    valid_labels = VALID_LABELS_LEGACY if legacy else VALID_LABELS_NEW

    auto = _auto_label(stl_paths, legacy)
    if auto is not None:
        click.echo("\n--- Patch labelling (auto from naming scheme) ---")
        for path, label in auto.items():
            click.echo(f"  {path.name}  ->  {label}")
        _echo_confirmed(auto, legacy)
        return auto

    if legacy:
        click.echo(
            "\n--- Patch labelling (legacy mode) ---\n"
            "Each file below is either the vessel wall surface or a capped opening.\n"
            "VMTK numbers caps geometrically; only you know which opening is the inlet.\n"
        )
    else:
        click.echo(
            "\n--- Patch labelling ---\n"
            "Each file below is the aneurysm sac surface, the parent vessel surface,\n"
            "or a capped opening (inlet or outlet).\n"
            "Labels: aneurysm_sac / parent_vessel / inlet / outlet\n"
        )

    labels: dict[Path, str] = {}
    for path in stl_paths:
        click.echo(f"File : {path.name}")
        click.echo(_describe(path))
        label = click.prompt(
            "  Label",
            type=click.Choice(valid_labels, case_sensitive=False),
        ).lower()
        labels[path] = label
        click.echo()

    errors = _validate_labels(labels, legacy)
    if errors:
        for e in errors:
            click.echo(f"ERROR: {e}", err=True)
        sys.exit(1)

    _echo_confirmed(labels, legacy)
    return labels


def _echo_confirmed(labels: dict[Path, str], legacy: bool) -> None:
    """Print the confirmed label assignment."""
    inlets = [p for p, l in labels.items() if l == "inlet"]
    outlets = [p for p, l in labels.items() if l == "outlet"]
    if legacy:
        walls = [p for p, l in labels.items() if l == "wall"]
        click.echo(
            f"Labels confirmed — wall: {walls[0].name}  "
            f"inlet: {inlets[0].name}  "
            f"outlets: {[p.name for p in outlets]}"
        )
    else:
        sacs = [p for p, l in labels.items() if l == "aneurysm_sac"]
        vessels = [p for p, l in labels.items() if l == "parent_vessel"]
        click.echo(
            f"Labels confirmed — aneurysm_sac: {sacs[0].name}  "
            f"parent_vessel: {vessels[0].name}  "
            f"inlet: {inlets[0].name}  "
            f"outlets: {[p.name for p in outlets]}"
        )
