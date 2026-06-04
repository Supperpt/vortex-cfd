"""
Phase B — Womersley pulsatile inlet profile.

Computes the exact Womersley analytical velocity profile for pulsatile flow in a
circular tube and writes it as OpenFOAM ``timeVaryingMappedFixedValue`` boundary data.

The Womersley solution decomposes Q(t) into Fourier harmonics and uses complex Bessel
functions to compute the frequency-dependent radial profile for each harmonic:

    U(r, t) = Σ_k Re[ C_k · W_k(r/R) · exp(i·k·ω·t) ] / A

where:
    k = 0  → Poiseuille (mean):  W_0 = 2·(1 - (r/R)²)
    k ≥ 1  → W_k(s) = 1 - J₀(Λ_k·s) / J₀(Λ_k),  s = r/R
    Λ_k = i^(3/2) · α_k,  α_k = R·√(k·ω/ν)  (Womersley number for harmonic k)
    C_k   = complex Fourier coefficient of Q(t) for harmonic k

scipy.special.jv(0, z) handles complex z directly — no polynomial approximations needed.

Approach confirmed by:
    https://github.com/JieWangnk/inlet-mapping-toolkit/blob/main/inlet_mapper/profiles.py

Delivery: ``timeVaryingMappedFixedValue``. Python writes per-face velocity vectors to
``constant/boundaryData/<patch>/`` after snappyHexMesh (face positions only exist then).
OpenFOAM reads and time-interpolates. No C++ compilation required.

Default inlet is always ``flowRateInletVelocity`` (parabolic). Womersley is opt-in only
via ``--womersley``.
"""

from pathlib import Path

import numpy as np
from scipy.special import jv

from .waveform import T_CYCLE


# ---------------------------------------------------------------------------
# Pure math — no OpenFOAM, no pyvista
# ---------------------------------------------------------------------------

def fourier_coefficients(waveform: np.ndarray, N_harm: int = 8) -> np.ndarray:
    """
    Compute the complex Fourier coefficients of the normalised flow-rate waveform.

    waveform : (N, 2) array of (t_norm, flow_norm) with flow_norm mean = 1.
    N_harm   : number of harmonics to retain (k = 1 … N_harm).

    Returns (N_harm + 1,) complex array where index k is the coefficient for
    harmonic k (k=0 is the DC / mean component).
    """
    q = waveform[:, 1]
    coeffs = np.fft.fft(q) / len(q)
    return coeffs[: N_harm + 1]


def _womersley_shape(r_norm: np.ndarray, alpha_k: float) -> np.ndarray:
    """
    Normalised Womersley shape function G_k(r/R) for Womersley number alpha_k.

        G_k(s) = W_k(s) / N_k

    where:
        W_k(s) = 1 - J₀(Λ·s) / J₀(Λ)      (unnormalised shape)
        N_k    = 1 - 2·J₁(Λ) / (Λ·J₀(Λ))   (normalisation so area-mean = 1)
        Λ      = i^(3/2) · alpha_k

    The area-weighted mean of Re(G_k) is 1, matching the Poiseuille G_0 = 2(1-s²).

    Returns complex (N_r,) array.
    """
    Lambda = (1j ** 1.5) * alpha_k
    J0_L = jv(0, Lambda)
    if np.abs(J0_L) < 1e-12:
        return np.ones_like(r_norm, dtype=complex)

    W_k = 1.0 - jv(0, Lambda * r_norm) / J0_L

    # Normalisation factor: area-weighted mean of W_k = N_k
    # 2 ∫₀¹ W_k(s) s ds = 1 - 2·J₁(Λ)/(Λ·J₀(Λ))
    N_k = 1.0 - 2.0 * jv(1, Lambda) / (Lambda * J0_L)
    if np.abs(N_k) < 1e-12:
        return np.ones_like(r_norm, dtype=complex)

    return W_k / N_k


def womersley_velocities(
    face_centers: np.ndarray,
    centroid: np.ndarray,
    normal: np.ndarray,
    radius: float,
    mean_velocity: float,
    waveform: np.ndarray,
    nu: float = 3.3e-6,
    t_cycle: float = T_CYCLE,
    N_harm: int = 8,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute the Womersley velocity vector at every face centre for each time step.

    face_centers : (N_faces, 3) face centroid positions [m]
    centroid     : (3,) inlet patch centroid [m]
    normal       : (3,) unit inward normal (pointing into the domain)
    radius       : float, inlet radius [m]
    mean_velocity: float, time-averaged inlet velocity [m/s]
    waveform     : (N_t, 2) normalised waveform (t_norm, flow_norm), mean = 1
    nu           : kinematic viscosity [m²/s]
    t_cycle      : cardiac period [s]
    N_harm       : number of Fourier harmonics

    Returns (t_abs, U_vectors):
        t_abs     : (N_t,) absolute times [s]
        U_vectors : (N_t, N_faces, 3) velocity vectors [m/s]
    """
    face_centers = np.asarray(face_centers)
    centroid = np.asarray(centroid)
    normal = np.asarray(normal, dtype=float)
    normal = normal / np.linalg.norm(normal)

    N_t = len(waveform)
    t_abs = waveform[:, 0] * t_cycle                          # (N_t,)
    omega = 2.0 * np.pi / t_cycle

    # Radial distance of each face from the inlet centroid (projected in the
    # plane of the inlet, i.e. perpendicular to the inward normal).
    delta = face_centers - centroid                            # (N_faces, 3)
    r = np.sqrt(np.sum(delta ** 2, axis=1)
                - np.dot(delta, normal) ** 2)                 # (N_faces,)
    r = np.clip(r, 0.0, radius)
    r_norm = r / radius                                        # (N_faces,)

    # Fourier coefficients of Q(t) scaled by Q_mean = mean_velocity * A
    # (area cancels in the normalised waveform, handled below per-harmonic)
    coeffs = fourier_coefficients(waveform, N_harm)            # (N_harm+1,)

    # Build the velocity magnitude profile U_mag[i_t, i_face]
    U_mag = np.zeros((N_t, len(face_centers)), dtype=float)

    # k = 0: steady Poiseuille component
    # U_poiseuille(r) = 2·U_mean·(1 - (r/R)²) — area-weighted mean = U_mean
    poiseuille = 2.0 * mean_velocity * (1.0 - r_norm ** 2)    # (N_faces,)
    U_mag += np.real(coeffs[0]) * poiseuille[np.newaxis, :]   # DC is real

    # k ≥ 1: oscillatory harmonics
    for k in range(1, N_harm + 1):
        if k >= len(coeffs):
            break
        alpha_k = radius * np.sqrt(k * omega / nu)
        W_k = _womersley_shape(r_norm, alpha_k)               # (N_faces,) complex
        phase = np.exp(1j * k * omega * t_abs)                 # (N_t,) complex
        # C_k * W_k(r) * exp(i k ω t) — C_k already scaled so mean flow = C_k * A
        # We work with velocities directly: C_k [m³/s] / A [m²] = C_k_vel [m/s]
        # But waveform is normalised (mean=1) and we already used mean_velocity in
        # the Poiseuille term; the harmonic contribution is:
        #   Re[ 2·C_k·W_k(r)·exp(ikωt) ]  (factor 2 for two-sided FFT)
        contrib = 2.0 * np.real(
            coeffs[k] * W_k[np.newaxis, :] * phase[:, np.newaxis]
        ) * mean_velocity                                      # (N_t, N_faces)
        U_mag += contrib

    # Brief retrograde flow near the wall is physically correct for Womersley
    # flow — do NOT clip. OpenFOAM handles both inflow and outflow at the same
    # patch boundary (the solver sees the net flux, and near-wall retrograde
    # flow is simply part of the pulsatile cycle).  Any clipping here would
    # inflate the computed area-mean above the waveform value at diastole.

    # Velocity vectors: U_mag × inward normal direction
    U_vectors = U_mag[:, :, np.newaxis] * normal[np.newaxis, np.newaxis, :]

    return t_abs, U_vectors


# ---------------------------------------------------------------------------
# OpenFOAM boundaryData file writer
# ---------------------------------------------------------------------------

_OF_HEADER = """\
FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    object      {obj};
}}
"""


def _write_foam_file(path: Path, cls: str, obj: str, body: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_OF_HEADER.format(cls=cls, obj=obj) + body, encoding="utf-8")


def write_boundary_data(
    case_dir: Path,
    patch_name: str,
    times: np.ndarray,
    face_centers: np.ndarray,
    U_vectors: np.ndarray,
) -> None:
    """
    Write the ``constant/boundaryData/<patch_name>/`` directory expected by
    ``timeVaryingMappedFixedValue``.

    Layout:
        constant/boundaryData/<patch>/
            points          ← face centroid positions (pointField)
            0.00000/U       ← velocity vectors at t=0
            0.01714/U       ← ... one directory per time step
            ...

    times      : (N_t,) absolute times [s]
    face_centers: (N_faces, 3) face centroid positions [m]
    U_vectors  : (N_t, N_faces, 3) velocity vectors [m/s]
    """
    case_dir = Path(case_dir)
    bd_dir = case_dir / "constant" / "boundaryData" / patch_name

    # Write points file
    n_pts = len(face_centers)
    pts_lines = "\n".join(
        f"( {p[0]:.10g} {p[1]:.10g} {p[2]:.10g} )" for p in face_centers
    )
    _write_foam_file(
        bd_dir / "points",
        cls="vectorField",
        obj="points",
        body=f"\n{n_pts}\n(\n{pts_lines}\n)\n",
    )

    # Write one U file per time step
    for i, t in enumerate(times):
        t_str = f"{t:.5f}"
        u_lines = "\n".join(
            f"( {v[0]:.10g} {v[1]:.10g} {v[2]:.10g} )" for v in U_vectors[i]
        )
        _write_foam_file(
            bd_dir / t_str / "U",
            cls="vectorField",
            obj="U",
            body=f"\n{n_pts}\n(\n{u_lines}\n)\n",
        )
