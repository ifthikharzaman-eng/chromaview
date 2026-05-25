"""Peak analysis: height extraction and signal-to-noise ratio.

Peak heights at called positions tell us about signal strength.
Signal-to-noise is computed as the ratio of the called base's peak height
to the mean of the other three channels at that position.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core.models import Chromatogram


@dataclass
class PeakInfo:
    """Per-base peak statistics."""

    position: int        # index into basecalls
    base: str
    peak_location: int   # sample point index
    peak_height: float   # called channel intensity at peak
    noise_height: float  # mean of other channels at peak
    snr: float           # signal-to-noise ratio


def peak_heights(chrom: Chromatogram) -> list[PeakInfo]:
    """Extract peak height and SNR for every called base."""
    results = []
    all_bases = ["A", "T", "G", "C"]

    for i, (base, loc) in enumerate(zip(chrom.basecalls, chrom.peak_locations)):
        loc = int(loc)
        base = base.upper()

        # Get signal at peak position for each channel
        heights = {}
        for b in all_bases:
            trace = chrom.traces.get(b)
            if trace is not None and loc < len(trace):
                heights[b] = float(trace[loc])
            else:
                heights[b] = 0.0

        signal = heights.get(base, 0.0)
        others = [heights[b] for b in all_bases if b != base]
        noise = np.mean(others) if others else 1.0
        snr = signal / noise if noise > 0 else float("inf")

        results.append(PeakInfo(
            position=i,
            base=base,
            peak_location=loc,
            peak_height=signal,
            noise_height=noise,
            snr=snr,
        ))
    return results


def signal_to_noise(chrom: Chromatogram) -> np.ndarray:
    """Compute per-base SNR as a 1-D array."""
    peaks = peak_heights(chrom)
    return np.array([p.snr for p in peaks])


def average_snr(chrom: Chromatogram) -> float:
    """Overall average signal-to-noise ratio for the trimmed region."""
    snr = signal_to_noise(chrom)
    start, end = chrom.trim_start, chrom.trim_end
    trimmed = snr[start:end]
    if len(trimmed) == 0:
        return 0.0
    # Use median to be robust against outlier peaks
    return float(np.median(trimmed))
