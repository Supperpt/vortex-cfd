"""
Womersley pulsatile inlet profile (opt-in via ``--womersley``).

The default inlet is a parabolic ``flowRateInletVelocity`` modulated by the
cardiac waveform.  That is fine for bulk flow, but the *shape* of a pulsatile
profile in a rigid tube is frequency-dependent: at high Womersley number the
core lags the pressure gradient and near-wall flow reverses before the core
does.  Publication-grade WSS needs the exact analytical solution.

The waveform is decomposed into Fourier harmonics; each harmonic k gets a
complex-Bessel radial shape:

    U(r, t) = U_mean · [ C_0 · 2(1 - s²)  +  2·Re( Σ_{k>=1} C_k · G_k(s) · e^{ikωt} ) ]

    s      = r/R
    G_k(s) = [1 - J0(Λs)/J0(Λ)] / N_k,   Λ = i^(3/2)·α_k,   α_k = R·sqrt(kω/ν)
    N_k    = 1 - 2·J1(Λ)/(Λ·J0(Λ))       (area-normalisation, so mean(G_k) = 1)

``scipy.special.jv`` takes complex arguments directly, so no polynomial
approximation is needed.  Each harmonic is area-normalised, which makes the
area-mean of U(r,t) reproduce ``U_mean · waveform(t)`` exactly — that identity
is what the unit tests assert.

**Delivery.**  ``timeVaryingMappedFixedValue`` reading per-face vectors from
``constant/boundaryData/inlet/``.  The alternative, ``codedFixedValue``, was
rejected: it needs runtime C++ compilation (fragile across ESI versions and on
locked-down HPC nodes), OpenFOAM has no native complex Bessel, and the maths
would stop being testable in Python.

**Multi-cycle.**  Unlike ``flowRateInletVelocity``'s ``outOfBounds repeat``,
``timeVaryingMappedFixedValue`` does *not* wrap around: past the last supplied
time it holds the final value forever.  Data must therefore be written for the
whole run, so ``womersley_velocities`` tiles the normalised cycle across all
``cycles`` and adds a terminal sample at exactly ``cycles·T``.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.special import jv

from .constants import NU
from .waveform import T_CYCLE

# Harmonics retained from the FFT.  8 resolves the ICA waveform's systolic
# upstroke; beyond ~10 the coefficients are numerical noise on a 100-point curve.
N_HARMONICS = 8


# ---------------------------------------------------------------------------
# Pure maths — no OpenFOAM, no pyvista
# ---------------------------------------------------------------------------

def fourier_coefficients(waveform: np.ndarray, n_harm: int = N_HARMONICS) -> np.ndarray:
    """
    Complex Fourier coefficients of the normalised flow-rate waveform.

    waveform : (N, 2) array of (t_norm, flow_norm), flow_norm mean = 1.
    Returns  : (n_harm + 1,) complex array; index 0 is the DC (mean) term.
    """
    q = np.asarray(waveform)[:, 1]
    coeffs = np.fft.fft(q) / len(q)
    return coeffs[: n_harm + 1]


def _womersley_shape(r_norm: np.ndarray, alpha_k: float) -> np.ndarray:
    """
    Area-normalised Womersley shape G_k(s) for Womersley number ``alpha_k``.

    Normalised so the area-weighted mean of G_k is 1, matching the Poiseuille
    term 2(1 - s²).  That is what lets the harmonics simply sum to the waveform.
    """
    lam = (1j ** 1.5) * alpha_k
    j0_lam = jv(0, lam)
    if np.abs(j0_lam) < 1e-12:
        # Λ sits on a zero of J0; the shape is degenerate, fall back to plug flow.
        return np.ones_like(r_norm, dtype=complex)

    shape = 1.0 - jv(0, lam * r_norm) / j0_lam

    # 2∫₀¹ shape(s)·s ds = 1 - 2·J1(Λ)/(Λ·J0(Λ))
    norm = 1.0 - 2.0 * jv(1, lam) / (lam * j0_lam)
    if np.abs(norm) < 1e-12:
        return np.ones_like(r_norm, dtype=complex)
    return shape / norm


def womersley_velocities(
    face_centers: np.ndarray,
    centroid: np.ndarray,
    normal: np.ndarray,
    radius: float,
    mean_velocity: float,
    waveform: np.ndarray,
    nu: float = NU,
    t_cycle: float = T_CYCLE,
    cycles: int = 1,
    n_harm: int = N_HARMONICS,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Velocity vector at every inlet face centre, for every time step of the run.

    normal : unit inward normal (into the domain).
    cycles : number of cardiac cycles to cover.  The normalised cycle is tiled
             and a terminal sample is appended at exactly ``cycles·t_cycle``,
             because timeVaryingMappedFixedValue does not repeat its data.

    Returns (times, U_vectors) with shapes (N_t,) and (N_t, N_faces, 3).
    """
    face_centers = np.asarray(face_centers, dtype=float)
    centroid = np.asarray(centroid, dtype=float)
    normal = np.asarray(normal, dtype=float)
    normal = normal / np.linalg.norm(normal)
    if radius <= 0:
        raise ValueError(f"inlet radius must be positive (got {radius})")
    if cycles < 1:
        raise ValueError(f"cycles must be >= 1 (got {cycles})")

    waveform = np.asarray(waveform, dtype=float)
    t_norm = waveform[:, 0]

    # Radial distance of each face from the axis, measured IN the inlet plane.
    delta = face_centers - centroid
    axial = delta @ normal
    r = np.linalg.norm(delta - np.outer(axial, normal), axis=1)
    r_norm = np.clip(r, 0.0, radius) / radius

    omega = 2.0 * np.pi / t_cycle
    coeffs = fourier_coefficients(waveform, n_harm)

    # Profile for one normalised cycle; every cycle repeats it, so evaluate once.
    n_t = len(waveform)
    u_mag = np.zeros((n_t, len(face_centers)), dtype=float)

    # k = 0 — steady Poiseuille, area-mean = mean_velocity.
    poiseuille = 2.0 * mean_velocity * (1.0 - r_norm ** 2)
    u_mag += float(np.real(coeffs[0])) * poiseuille[np.newaxis, :]

    # k >= 1 — oscillatory harmonics. Factor 2 accounts for the conjugate half
    # of the two-sided FFT.
    t_in_cycle = t_norm * t_cycle
    for k in range(1, min(n_harm, len(coeffs) - 1) + 1):
        alpha_k = radius * np.sqrt(k * omega / nu)
        shape = _womersley_shape(r_norm, alpha_k)
        phase = np.exp(1j * k * omega * t_in_cycle)
        u_mag += 2.0 * mean_velocity * np.real(
            coeffs[k] * shape[np.newaxis, :] * phase[:, np.newaxis]
        )

    # Near-wall retrograde flow during diastole is physically correct for
    # Womersley flow and must NOT be clipped: clipping would raise the
    # area-mean above the waveform value exactly where the waveform is lowest.

    # Tile across the run, then close the interval so the solver never runs past
    # the supplied data and silently freezes the inlet.
    times = np.concatenate(
        [t_in_cycle + c * t_cycle for c in range(cycles)] + [[cycles * t_cycle]]
    )
    profiles = np.concatenate([u_mag] * cycles + [u_mag[:1]], axis=0)

    vectors = profiles[:, :, np.newaxis] * normal[np.newaxis, np.newaxis, :]
    return times, vectors


# ---------------------------------------------------------------------------
# OpenFOAM boundaryData writer
# ---------------------------------------------------------------------------

_HEADER = """\
FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    object      {obj};
}}
"""


def _write_vector_field(path: Path, obj: str, vectors: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    body = "\n".join(f"({v[0]:.10g} {v[1]:.10g} {v[2]:.10g})" for v in vectors)
    path.write_text(
        _HEADER.format(cls="vectorField", obj=obj)
        + f"\n{len(vectors)}\n(\n{body}\n)\n",
        encoding="utf-8",
    )


def write_boundary_data(
    case_dir: Path,
    patch_name: str,
    times: np.ndarray,
    face_centers: np.ndarray,
    vectors: np.ndarray,
) -> Path:
    """
    Write ``constant/boundaryData/<patch_name>/`` for timeVaryingMappedFixedValue:

        points        face centre positions
        <t>/U         velocity vectors at each supplied time

    Returns the boundaryData directory.
    """
    case_dir = Path(case_dir)
    if len(times) != len(vectors):
        raise ValueError(
            f"times and vectors disagree: {len(times)} times, {len(vectors)} fields"
        )
    bd_dir = case_dir / "constant" / "boundaryData" / patch_name

    _write_vector_field(bd_dir / "points", "points", np.asarray(face_centers))
    for t, field in zip(times, vectors):
        _write_vector_field(bd_dir / f"{t:.5f}" / "U", "U", field)
    return bd_dir
