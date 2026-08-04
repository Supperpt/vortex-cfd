# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/), and this project adheres to
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- Neck inflow metrics (`--neck-metrics`, opt-in and unvalidated): neck inflow rate,
  net flux and peak velocity, computed by slicing the volume velocity field at the
  aneurysm neck orifice. The orifice is fitted to the open boundary loop of
  `aneurysm_sac.stl` and recorded in `neck_plane_resolved.json`.

### Changed
- The neck `surfaceFieldValue` function objects have been removed rather than
  re-enabled. The clinical neck inflow rate is the integral of the *inward* part of
  the velocity only, which no function-object operation can express;
  `areaNormalIntegrate` measures net flux, which averages to ~0 through a sealed sac.
- `surfaceFieldValue` parsing now strips OpenFOAM's parentheses from vector values,
  warns about unparsable rows instead of dropping them silently, and takes an explicit
  vector-reduction argument so a signed quantity cannot lose its sign by default.

### Fixed
- BUG-010: `neck_peak_velocity_ms` was always null — the parser could not read
  parenthesised vector tokens, and the failure was swallowed silently.
- BUG-011: `neck_peak_flow_rate_m3s` used a sign-naive `max()` on signed flux, which
  returned the *smallest*-magnitude value whenever flux was negative all cycle.
- CAVEAT-012: the neck sampling plane was infinite and integrated the whole
  parent-vessel cross-section instead of the sac orifice.

## [1.0.0] — 2026-06-17

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

[1.0.0]: https://github.com/Supperpt/vortex-cfd/releases/tag/v1.0.0
