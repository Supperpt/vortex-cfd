"""Default and user-supplied pulsatile cardiac waveforms."""

import csv
from pathlib import Path

import numpy as np

T_CYCLE = 0.857   # seconds — 70 bpm
N_POINTS = 100    # samples per cycle in the flow-rate table


def _default_shape(t_norm: np.ndarray) -> np.ndarray:
    """
    Normalised analytical ICA waveform (mean = 1, systolic peak at ~30 % of cycle).
    Built from 3 Fourier harmonics tuned to match the typical internal-carotid shape.
    Negative excursions (brief retrograde flow) are clipped to 2 % of the mean,
    then the array is renormalised so its mean is exactly 1.0.
    """
    w = 2 * np.pi
    shape = (
        1.0
        + 1.2 * np.cos(w * t_norm - 1.6)
        + 0.20 * np.cos(2 * w * t_norm - 0.7)
        + 0.10 * np.cos(3 * w * t_norm - 0.3)
    )
    shape = np.maximum(shape, 0.02)
    return shape / np.mean(shape)


def load_waveform(csv_path: str | None) -> np.ndarray:
    """
    Return an (N, 2) array of (t_normalised, flow_normalised) for one cycle.
    flow_normalised has mean = 1 over the interval.

    If csv_path is None the built-in analytical shape is used.
    The CSV must have two columns: time_normalised (0–1) and flow_normalised.
    Lines starting with '#' and blank lines are ignored.
    """
    t_norm = np.linspace(0, 1, N_POINTS, endpoint=False)

    if csv_path is None:
        return np.column_stack([t_norm, _default_shape(t_norm)])

    rows: list[tuple[float, float]] = []
    with open(csv_path, newline="") as fh:
        for row in csv.reader(fh):
            if not row or row[0].startswith("#"):
                continue
            try:
                rows.append((float(row[0]), float(row[1])))
            except (ValueError, IndexError):
                continue

    if len(rows) < 2:
        raise ValueError(f"Waveform CSV '{csv_path}' must contain at least 2 data rows.")

    t_arr = np.array([r[0] for r in rows])
    q_arr = np.array([r[1] for r in rows])
    q_arr = q_arr / np.mean(q_arr)
    return np.column_stack([t_arr, q_arr])
