"""Interactive wall / inlet / outlet labelling of VORTEX STL caps."""
from __future__ import annotations

import sys
from pathlib import Path

import click
import pyvista as pv

# Legacy mode: single wall patch.
VALID_LABELS_LEGACY = ("wall", "inlet", "outlet")

# New (default) mode: aneurysm sac + parent vessel as separate wall patches.
VALID_LABELS_NEW = ("aneurysm_sac", "parent_vessel", "inlet", "outlet")


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
    """
    valid_labels = VALID_LABELS_LEGACY if legacy else VALID_LABELS_NEW

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

    if errors:
        for e in errors:
            click.echo(f"ERROR: {e}", err=True)
        sys.exit(1)

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

    return labels
