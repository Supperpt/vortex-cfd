"""Click-based entry point for vortex-cfd."""

import sys
from pathlib import Path

import click

from .env_check import check_openfoam
from .patch_labeller import label_patches
from .scaling import scale_stls
from .waveform import load_waveform
from .case_builder import build_case
from .runner import run_pipeline


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--stl-dir",       required=True, type=click.Path(exists=True, file_okay=False),
              help="Directory containing wall + cap STLs from VORTEX (--split-patches).")
@click.option("--cycles",        required=True, type=int, default=3, show_default=True,
              help="Number of cardiac cycles to simulate (first discarded, last analysed).")
@click.option("--mean-velocity", required=True, type=float,
              help="Time-averaged inlet velocity in m/s (typical ICA: 0.3–0.5).")
@click.option("--waveform",      "waveform_csv", default=None,
              type=click.Path(exists=True, dir_okay=False),
              help="Optional 2-column CSV (time_norm, flow_norm) for the pulse shape.")
@click.option("--cores",         default=None, type=int,
              help="CPU cores for parallel meshing/solving (default: all available).")
@click.option("--out-dir",       default=".", show_default=True,
              type=click.Path(file_okay=False),
              help="Parent directory for the OpenFOAM case directory.")
def main(stl_dir, cycles, mean_velocity, waveform_csv, cores, out_dir):
    """
    Automated pulsatile CFD for cerebral aneurysms.

    Takes the split-patch STL output of VORTEX (wall + capped openings) and
    produces a complete, runnable OpenFOAM case with WSS-resolved boundary layers.
    """
    # 1. Validate OpenFOAM environment
    of_env = check_openfoam()
    click.echo(f"OpenFOAM {of_env['version']} detected at {of_env['root'] or '(sourced)'}")

    # 2. Find STLs
    stl_paths = sorted(Path(stl_dir).glob("*.stl"))
    if not stl_paths:
        click.echo(f"ERROR: No STL files found in {stl_dir}", err=True)
        sys.exit(1)
    click.echo(f"Found {len(stl_paths)} STL file(s): {[p.name for p in stl_paths]}")

    # 3. Interactive labelling
    labels = label_patches(stl_paths)

    # 4. Scale mm → m if needed and assign canonical names
    scaled_stls = scale_stls(stl_paths, labels)

    # 5. Load waveform
    waveform = load_waveform(waveform_csv)
    if waveform_csv:
        click.echo(f"Using user waveform: {waveform_csv}")
    else:
        click.echo("Using built-in analytical ICA waveform.")

    # 6. Build OpenFOAM case directory
    case_dir = build_case(
        scaled_stls=scaled_stls,
        labels=labels,
        cycles=cycles,
        mean_velocity=mean_velocity,
        waveform=waveform,
        cores=cores,
        out_dir=out_dir,
    )
    click.echo(f"Case directory created: {case_dir}")

    # 7. Run the pipeline
    run_pipeline(case_dir=case_dir, of_env=of_env, cores=cores)
