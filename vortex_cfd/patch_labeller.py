"""Interactive wall / inlet / outlet labelling of VORTEX STL caps."""

import sys
from pathlib import Path

import click
import pyvista as pv

VALID_LABELS = ("wall", "inlet", "outlet")


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


def label_patches(stl_paths: list[Path]) -> dict[Path, str]:
    """
    Prompt the user to assign wall / inlet / outlet to each STL.
    Validates that the result has exactly one wall and one inlet.
    Returns a {Path: label} dict.
    """
    click.echo(
        "\n--- Patch labelling ---\n"
        "Each file below is either the vessel wall surface or a capped opening.\n"
        "VMTK numbers caps geometrically; only you know which opening is the inlet.\n"
    )

    labels: dict[Path, str] = {}
    for path in stl_paths:
        click.echo(f"File : {path.name}")
        click.echo(_describe(path))
        label = click.prompt(
            "  Label",
            type=click.Choice(VALID_LABELS, case_sensitive=False),
        ).lower()
        labels[path] = label
        click.echo()

    walls = [p for p, l in labels.items() if l == "wall"]
    inlets = [p for p, l in labels.items() if l == "inlet"]
    outlets = [p for p, l in labels.items() if l == "outlet"]

    errors: list[str] = []
    if len(walls) != 1:
        errors.append(f"Expected exactly 1 wall, got {len(walls)}: {[p.name for p in walls]}")
    if len(inlets) != 1:
        errors.append(f"Expected exactly 1 inlet, got {len(inlets)}: {[p.name for p in inlets]}")
    if not outlets:
        errors.append("Expected at least 1 outlet.")

    if errors:
        for e in errors:
            click.echo(f"ERROR: {e}", err=True)
        sys.exit(1)

    click.echo(
        f"Labels confirmed — wall: {walls[0].name}  "
        f"inlet: {inlets[0].name}  "
        f"outlets: {[p.name for p in outlets]}"
    )
    return labels
