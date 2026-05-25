"""Quality analysis and trimming for chromatogram data.

Implements the modified Mott trimming algorithm used by phred/Chromas:
    1. Subtract a quality threshold from each Phred score.
    2. Compute a running sum of these adjusted values.
    3. Find the region where the running sum is maximized — that region
       represents the longest stretch of above-threshold quality.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core.models import Chromatogram


def trim_low_quality(
    chrom: Chromatogram,
    threshold: int = 20,
    window: int = 10,
    apply: bool = True,
) -> tuple[int, int]:
    """Determine trim boundaries using the modified Mott algorithm.

    Args:
        chrom: Chromatogram to analyze.
        threshold: Phred score cutoff (default 20 = 1% error rate).
        window: Minimum number of consecutive good bases at edges.
        apply: If True, set chrom.trim_start and chrom.trim_end.

    Returns:
        (trim_start, trim_end) indices into basecalls.
    """
    scores = chrom.quality_scores
    if len(scores) == 0:
        return (0, 0)

    n = len(scores)

    # Modified Mott: subtract threshold, then find max-sum subarray
    adjusted = scores.astype(float) - threshold
    cumsum = np.cumsum(adjusted)

    # Find the subarray [i..j] with maximum sum (Kadane's variant)
    max_sum = float("-inf")
    best_start, best_end = 0, n
    current_start = 0
    running = 0.0

    for i in range(n):
        running += adjusted[i]
        if running > max_sum:
            max_sum = running
            best_start = current_start
            best_end = i + 1
        if running < 0:
            running = 0.0
            current_start = i + 1

    # Additional edge refinement: skip leading/trailing bases below threshold
    while best_start < best_end and scores[best_start] < threshold:
        best_start += 1
    while best_end > best_start and scores[best_end - 1] < threshold:
        best_end -= 1

    # Ensure minimum window of good bases at boundaries
    if best_end - best_start < window:
        best_start = 0
        best_end = n

    if apply:
        chrom.trim_start = best_start
        chrom.trim_end = best_end

    return (best_start, best_end)


@dataclass
class QualitySummary:
    """Aggregate quality statistics for a chromatogram."""

    total_bases: int
    trimmed_bases: int
    mean_quality: float
    median_quality: float
    q20_count: int      # bases with Phred >= 20
    q30_count: int      # bases with Phred >= 30
    q20_fraction: float
    q30_fraction: float
    min_quality: int
    max_quality: int


def quality_summary(chrom: Chromatogram) -> QualitySummary:
    """Compute aggregate quality statistics."""
    qs = chrom.trimmed_quality
    n = len(qs)
    if n == 0:
        return QualitySummary(
            total_bases=chrom.num_bases,
            trimmed_bases=0,
            mean_quality=0.0,
            median_quality=0.0,
            q20_count=0, q30_count=0,
            q20_fraction=0.0, q30_fraction=0.0,
            min_quality=0, max_quality=0,
        )

    return QualitySummary(
        total_bases=chrom.num_bases,
        trimmed_bases=n,
        mean_quality=float(np.mean(qs)),
        median_quality=float(np.median(qs)),
        q20_count=int(np.sum(qs >= 20)),
        q30_count=int(np.sum(qs >= 30)),
        q20_fraction=float(np.sum(qs >= 20)) / n,
        q30_fraction=float(np.sum(qs >= 30)) / n,
        min_quality=int(np.min(qs)),
        max_quality=int(np.max(qs)),
    )
