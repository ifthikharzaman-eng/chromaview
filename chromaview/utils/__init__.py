"""Shared utility functions."""
from __future__ import annotations

from pathlib import Path


def detect_file_format(filepath: str | Path) -> str:
    """Detect chromatogram file format from extension."""
    ext = Path(filepath).suffix.lower()
    if ext == ".ab1":
        return "ab1"
    elif ext == ".scf":
        return "scf"
    else:
        raise ValueError(f"Unsupported file extension: {ext}")


def parse_fasta_string(text: str) -> str:
    """Extract raw sequence from a FASTA-formatted string."""
    lines = text.strip().split("\n")
    seq_lines = [l.strip() for l in lines if not l.startswith(">") and l.strip()]
    return "".join(seq_lines).replace(" ", "").upper()
