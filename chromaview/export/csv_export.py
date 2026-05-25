"""CSV export for base-by-base data."""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Union

from ..core.models import Chromatogram


def export_csv(
    chrom: Chromatogram,
    filepath: Union[str, Path],
    trimmed: bool = True,
    include_traces: bool = False,
) -> Path:
    """Export per-base data to CSV.

    Columns: Position, Base, Quality, PeakLocation[, A_height, T_height, G_height, C_height]

    Args:
        chrom: Source chromatogram.
        filepath: Output path.
        trimmed: Export only the trimmed region.
        include_traces: Add trace heights at each peak position.

    Returns:
        Path to the written file.
    """
    filepath = Path(filepath)
    start = chrom.trim_start if trimmed else 0
    end = chrom.trim_end if trimmed else chrom.num_bases

    headers = ["Position", "Base", "Quality", "PeakLocation"]
    if include_traces:
        headers.extend(["A_height", "T_height", "G_height", "C_height"])

    with open(filepath, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)

        for i in range(start, end):
            loc = int(chrom.peak_locations[i]) if i < len(chrom.peak_locations) else 0
            quality = int(chrom.quality_scores[i]) if i < len(chrom.quality_scores) else 0
            row = [i + 1, chrom.basecalls[i], quality, loc]

            if include_traces:
                for base in ["A", "T", "G", "C"]:
                    trace = chrom.traces.get(base)
                    if trace is not None and loc < len(trace):
                        row.append(f"{trace[loc]:.1f}")
                    else:
                        row.append("0.0")

            writer.writerow(row)

    return filepath
