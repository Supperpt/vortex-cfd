"""Generate the bundled default ICA inlet waveform: ``ica_ford2005.csv``.

This is a one-time / reproducible data-generation script — it is NOT imported at
runtime (the pipeline ships and reads the static CSV it produces).  Run it only
to regenerate the CSV if the digitisation ever needs revisiting.

Source
------
Ford MD, Alperin N, Lee SH, Holdsworth DW, Steinman DA.
"Characterization of volumetric flow rate waveforms in the normal internal
carotid and vertebral arteries." Physiol Meas 2005;26(4):477-488.
DOI 10.1088/0967-3334/26/4/013.

The archetypal ICA waveform feature points are taken verbatim from the paper's
**Table 2** (ICA columns).  Timings are in ms relative to feature point H0;
amplitudes are already normalised to the cycle-averaged flow rate (VFR/VFR_avg),
so the cycle mean is ~1 by construction.  M0 (the global minimum) is the natural
cycle start; the mean R-R interval reported by the paper is ~885 ms, so the next
beat's M0 falls at ~885 ms after this beat's M0 and closes the cycle.

We fit a *periodic* cubic spline through the feature points (the paper itself fit
the feature points with a cubic spline), resample onto a uniform 100-point grid
over t_norm in [0, 1) (endpoint excluded so OpenFOAM's ``outOfBounds repeat`` is
seamless), and renormalise so the cycle mean is exactly 1.0.
"""

from pathlib import Path

import numpy as np

# --- Ford et al. (2005) Table 2, ICA archetypal waveform feature points -------
# (feature, timing_ms_relative_to_H0, amplitude_VFR_over_VFRavg)
FEATURE_POINTS = [
    ("M0", -61, 0.68),
    ("H0",   0, 1.18),
    ("P1",  45, 1.66),
    ("H1",  84, 1.42),
    ("M1", 141, 1.20),
    ("P2", 187, 1.22),
    ("H2", 241, 1.08),
    ("M2", 280, 0.94),
    ("H3", 312, 1.02),
    ("P3", 350, 1.09),
    ("D1", 468, 0.94),
    ("D2", 586, 0.83),
    ("D3", 704, 0.76),
    ("D4", 822, 0.68),
]
RR_INTERVAL_MS = 885.0   # mean R-R interval reported by Ford et al. (68 bpm)
N_POINTS = 100           # must match vortex_cfd.waveform.N_POINTS


def _periodic_cubic_spline(x, y, period, x_eval):
    """Evaluate a periodic natural cubic spline through (x, y) at x_eval.

    x must be strictly increasing within one period; the curve wraps so that
    x[0]+period maps back to y[0].  Pure-numpy (no scipy) solve of the cyclic
    second-derivative system.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    n = len(x)
    # interval lengths, with the last interval wrapping across the period
    h = np.empty(n)
    h[:-1] = np.diff(x)
    h[-1] = (x[0] + period) - x[-1]

    # cyclic tridiagonal system  A M = b  for second derivatives M (M[n] == M[0])
    A = np.zeros((n, n))
    b = np.zeros(n)
    for i in range(n):
        im1, ip1 = (i - 1) % n, (i + 1) % n
        hm = h[im1]  # length of interval ending at knot i
        hi = h[i]    # length of interval starting at knot i
        A[i, im1] += hm
        A[i, i] += 2.0 * (hm + hi)
        A[i, ip1] += hi
        b[i] = 6.0 * ((y[ip1] - y[i]) / hi - (y[i] - y[im1]) / hm)
    M = np.linalg.solve(A, b)
    M = np.append(M, M[0])  # close the loop: M[n] == M[0]

    # locate each eval point in its interval (within [x[0], x[0]+period))
    xr = (x_eval - x[0]) % period + x[0]
    idx = np.searchsorted(x, xr, side="right") - 1
    idx = np.clip(idx, 0, n - 1)

    out = np.empty_like(xr)
    for k, (xx, i) in enumerate(zip(xr, idx)):
        hi = h[i]
        xi = x[i]
        xi1 = x[i] + hi  # next knot position (may exceed x[-1] on the wrap)
        a = (xi1 - xx) / hi
        c = (xx - xi) / hi
        yi1 = y[(i + 1) % n]
        out[k] = (
            M[i] * a**3 * hi**2 / 6.0
            + M[i + 1] * c**3 * hi**2 / 6.0
            + (y[i] - M[i] * hi**2 / 6.0) * a
            + (yi1 - M[i + 1] * hi**2 / 6.0) * c
        )
    return out


def build_waveform():
    """Return (t_norm, flow_norm) arrays of length N_POINTS, mean(flow_norm)=1."""
    # shift so the cycle starts at M0 = 0 ms; period is the M0->next-M0 interval
    t0 = FEATURE_POINTS[0][1]            # -61 ms
    x = np.array([t for _, t, _ in FEATURE_POINTS], float) - t0   # 0 .. 883 ms
    y = np.array([a for _, _, a in FEATURE_POINTS], float)
    period = RR_INTERVAL_MS              # 885 ms; next M0 closes the cycle

    t_norm = np.linspace(0.0, 1.0, N_POINTS, endpoint=False)
    flow = _periodic_cubic_spline(x, y, period, t_norm * period)
    flow = flow / flow.mean()           # enforce exact cycle mean = 1
    return t_norm, flow


def main():
    t_norm, flow = build_waveform()
    out = Path(__file__).with_name("ica_ford2005.csv")
    header = (
        "# Default ICA inlet waveform for vortex-cfd.\n"
        "# Ford MD, Alperin N, Lee SH, Holdsworth DW, Steinman DA. Characterization of\n"
        "# volumetric flow rate waveforms in the normal internal carotid and vertebral\n"
        "# arteries. Physiol Meas 2005;26(4):477-488. DOI 10.1088/0967-3334/26/4/013.\n"
        "# Periodic cubic spline through the Table 2 ICA archetypal feature points,\n"
        "# resampled to 100 points over t_norm in [0,1), renormalised to cycle mean = 1.\n"
        "# Columns: t_norm (0-1, cycle fraction), flow_norm (VFR/VFR_avg, mean = 1).\n"
        "# Regenerate with: python -m vortex_cfd.data.generate_ica_ford2005\n"
    )
    lines = [f"{t:.6f},{q:.6f}" for t, q in zip(t_norm, flow)]
    out.write_text(header + "\n".join(lines) + "\n")

    # diagnostics
    peak_i = int(np.argmax(flow))
    print(f"wrote {out}  ({len(flow)} rows)")
    print(f"  mean      = {flow.mean():.6f}  (target 1.0)")
    print(f"  min       = {flow.min():.4f} at t_norm={t_norm[flow.argmin()]:.3f}")
    print(f"  peak      = {flow.max():.4f} at t_norm={t_norm[peak_i]:.3f}")
    print(f"  peak/mean = {flow.max():.3f}  (Ford P1 = 1.66)")
    print(f"  all positive (antegrade): {bool(np.all(flow > 0))}")


if __name__ == "__main__":
    main()
