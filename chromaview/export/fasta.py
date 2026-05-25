"""FASTA sequence export."""
from __future__ import annotations

from pathlib import Path
from typing import Union

from ..core.models import Chromatogram


def export_fasta(
    chrom: Chromatogram,
    filepath: Union[str, Path],
    trimmed: bool = True,
    line_width: int = 80,
) -> Path:
    """Export the called sequence to FASTA format.

    Args:
        chrom: Source chromatogram.
        filepath: Output file path.
        trimmed: If True, export only the trimmed region.
        line_width: Characters per line in the sequence block.

    Returns:
        Path to the written file.
    """
    filepath = Path(filepath)
    name = chrom.metadata.sample_name or filepath.stem
    seq = chrom.trimmed_sequence if trimmed else chrom.sequence

    # Build header with useful metadata
    header_parts = [name]
    if trimmed:
        header_parts.append(f"trimmed:{chrom.trim_start+1}-{chrom.trim_end}")
    header_parts.append(f"length:{len(seq)}")
    if chrom.is_edited:
        header_parts.append(f"edited:{len(chrom.get_edits())}changes")
    header = " ".join(header_parts)

    # Wrap sequence at line_width
    lines = [f">{header}"]
    for i in range(0, len(seq), line_width):
        lines.append(seq[i : i + line_width])

    filepath.write_text("\n".join(lines) + "\n")
    return filepath
