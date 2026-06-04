"""Click-based entry point for vortex-cfd."""

import logging
import sys
from pathlib import Path

import click

log = logging.getLogger("vortex_cfd")

from .env_check import check_openfoam
from .patch_labeller import label_patches
from .scaling import scale_stls
from .waveform import load_waveform
from .case_builder import build_case
from .runner import run_pipeline, run_postprocess_only


def _add_file_handler(case_dir: Path) -> None:
    """Attach a FileHandler writing to <case_dir>/vortex_cfd.log."""
    fh = logging.FileHandler(case_dir / "vortex_cfd.log")
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    logging.getLogger("vortex_cfd").addHandler(fh)


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--stl-dir",       default=None, type=click.Path(exists=True, file_okay=False),
              help="Directory containing wall + cap STLs from VORTEX (--split-patches).")
@click.option("--cycles",        type=int, default=3, show_default=True,
              help="Number of cardiac cycles to simulate (first discarded, last analysed).")
@click.option("--mean-velocity", default=None, type=float,
              help="Time-averaged inlet velocity in m/s (typical ICA: 0.3–0.5).")
@click.option("--waveform",      "waveform_csv", default=None,
              type=click.Path(exists=True, dir_okay=False),
              help="Optional 2-column CSV (time_norm, flow_norm) for the pulse shape.")
@click.option("--cores",         default=None, type=int,
              help="CPU cores for parallel meshing/solving (default: all available).")
@click.option("--out-dir",       default=".", show_default=True,
              type=click.Path(file_okay=False),
              help="Parent directory for the OpenFOAM case directory.")
@click.option("--postprocess",   is_flag=True, default=False,
              help="Compute WSS/TAWSS/OSI biomarkers and write metrics_report.json "
                   "(adds the wallShearStress + fieldAverage function objects).")
@click.option("--postprocess-only", "postprocess_only", default=None,
              type=click.Path(exists=True, file_okay=False),
              help="Skip meshing/solving; compute biomarkers on an existing solved "
                   "case directory (re-uses or regenerates the WSS field).")
@click.option("--womersley",     is_flag=True, default=False,
              help="Use the full Womersley analytical inlet profile instead of "
                   "the default parabolic (flowRateInletVelocity). Requires a "
                   "second meshing step to write per-face velocity data. "
                   "Recommended for academic publication.")
@click.option("--dry-run",       is_flag=True, default=False,
              help="Build the case directory but skip meshing and solving. "
                   "Useful for inspecting generated dicts before committing CPU time.")
def main(stl_dir, cycles, mean_velocity, waveform_csv, cores, out_dir,
         postprocess, postprocess_only, womersley, dry_run):
    """
    Automated pulsatile CFD for cerebral aneurysms.

    Takes the split-patch STL output of VORTEX (wall + capped openings) and
    produces a complete, runnable OpenFOAM case with WSS-resolved boundary layers.
    """
    # Standalone post-processing of an existing case — no meshing/solving.
    if postprocess_only:
        of_env = check_openfoam()
        log.info("OpenFOAM %s detected at %s", of_env['version'], of_env['root'] or '(sourced)')
        log.info("Post-processing existing case: %s", postprocess_only)
        _add_file_handler(Path(postprocess_only))
        run_postprocess_only(Path(postprocess_only), of_env, cycles=cycles)
        return

    # Normal run requires the patient-specific inputs.
    if stl_dir is None or mean_velocity is None:
        log.error("--stl-dir and --mean-velocity are required "
                  "(unless using --postprocess-only).")
        sys.exit(1)

    # 1. Validate OpenFOAM environment
    of_env = check_openfoam()
    log.info("OpenFOAM %s detected at %s", of_env['version'], of_env['root'] or '(sourced)')

    # 2. Find STLs
    stl_paths = sorted(Path(stl_dir).glob("*.stl"))
    if not stl_paths:
        log.error("No STL files found in %s", stl_dir)
        sys.exit(1)
    log.info("Found %d STL file(s): %s", len(stl_paths), [p.name for p in stl_paths])

    # 3. Interactive labelling
    labels = label_patches(stl_paths)

    # 4. Scale mm → m if needed and assign canonical names
    scaled_stls = scale_stls(stl_paths, labels)

    # 5. Load waveform
    waveform = load_waveform(waveform_csv)
    if waveform_csv:
        log.info("Using user waveform: %s", waveform_csv)
    else:
        log.info("Using built-in analytical ICA waveform.")

    # 6. Build OpenFOAM case directory
    case_dir, inlet_params = build_case(
        scaled_stls=scaled_stls,
        labels=labels,
        cycles=cycles,
        mean_velocity=mean_velocity,
        waveform=waveform,
        cores=cores,
        out_dir=out_dir,
        postprocess=postprocess,
        womersley=womersley,
    )
    log.info("Case directory created: %s", case_dir)
    _add_file_handler(case_dir)

    if dry_run:
        log.info("--dry-run: case directory built; meshing and solving skipped.")
        return

    # 7. Run the pipeline
    run_pipeline(
        case_dir=case_dir,
        of_env=of_env,
        cores=cores,
        cycles=cycles,
        postprocess_metrics=postprocess,
        womersley=womersley,
        inlet_params=inlet_params,
    )
