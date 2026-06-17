# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/), and this project adheres to
[Semantic Versioning](https://semver.org/).

## [0.1.0] — 2026-06-17

First public release, accompanying the ARTERY26 proof-of-concept study.

### Added
- Single-command CT-to-CFD pipeline (`vortex-cfd`): semantic STL patch labelling,
  mm→m scaling, OpenFOAM case generation from Jinja2 templates, meshing
  (snappyHexMesh with prismatic boundary layers), and pulsatile `pimpleFoam` solve.
- Hemodynamic biomarker post-processing (`--postprocess` / `--postprocess-only`):
  TAWSS, OSI, low-WSS and high-OSI area fractions, with a `metrics_report.json` output.
- Aneurysm (two-patch) mode: biomarkers scoped to the aneurysm dome, with
  normalised WSS (sac/parent) and sac pressure (mean/peak).
- Default literature-standard Ford et al. (2005) ICA inflow waveform, with
  `--waveform` override.
- Auto-labelling of STL files from a default naming scheme for unattended batch runs.

### Known limitations
- **Neck-flow metrics are disabled** (neck inflow rate and peak velocity). The
  underlying function objects are commented out pending validation: the infinite
  sampling plane integrates the whole parent-vessel cross-section rather than the
  sac orifice, and the flux/peak-velocity parsing is unreliable. Re-enabling these
  is the next planned milestone. Validated outputs are TAWSS, OSI, normalised WSS,
  and sac pressure.

[0.1.0]: https://github.com/Supperpt/vortex-cfd/releases/tag/v0.1.0
