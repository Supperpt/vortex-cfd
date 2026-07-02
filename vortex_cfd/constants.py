"""Shared physical constants for the vortex-cfd pipeline.

Single source of truth for the fluid properties used both when building the
OpenFOAM case (``case_builder``) and when post-processing it (``postprocess``),
so the two can never drift out of sync.
"""
from __future__ import annotations

# Default fluid properties — whole blood.
RHO = 1060.0      # kg/m^3, density
NU = 3.3e-6       # m^2/s, kinematic viscosity
