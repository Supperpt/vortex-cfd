"""Default and user-supplied pulsatile cardiac waveforms."""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

T_CYCLE = 0.857   # seconds — 70 bpm
N_POINTS = 100    # samples per cycle in the flow-rate table

# Built-in default inlet waveform: the Ford et al. (2005) archetypal internal
# carotid artery (ICA) volumetric flow-rate shape, normalised to cycle mean = 1.
# Digitised from Table 2 of the paper — see data/generate_ica_ford2005.py for
# provenance and regeneration.
DEFAULT_WAVEFORM_CSV = Path(__file__).parent / "data" / "ica_ford2005.csv"


def _read_waveform_csv(csv_path) -> tuple[np.ndarray, np.ndarray]:
    """
    Parse a 2-column (t_normalised, flow_normalised) CSV.
    Lines starting with '#' and blank lines are ignored; unparseable rows skipped.
    """
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
    return t_arr, q_arr


def load_waveform(csv_path: str | None) -> np.ndarray:
    """
    Return an (N, 2) array of (t_normalised, flow_normalised) for one cycle.
    flow_normalised has mean = 1 over the interval.

    If csv_path is None the built-in Ford et al. (2005) ICA waveform is used
    (bundled at data/ica_ford2005.csv). A user CSV must have two columns:
    time_normalised (0–1) and flow_normalised. Lines starting with '#' and blank
    lines are ignored. Either way the flow column is renormalised to mean = 1, so
    the magnitude is set solely by --mean-velocity downstream.
    """
    path = DEFAULT_WAVEFORM_CSV if csv_path is None else csv_path
    t_arr, q_arr = _read_waveform_csv(path)
    q_arr = q_arr / np.mean(q_arr)
    return np.column_stack([t_arr, q_arr])
