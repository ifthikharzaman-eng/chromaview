"""Data models for chromatogram representation.

The Chromatogram class is the central data structure. It holds raw trace data,
base calls, quality scores, and metadata. All modules operate on this model,
keeping parsing, visualization, and analysis cleanly separated.

Design decisions:
- NumPy arrays for trace data (fast slicing, compatible with PyQtGraph)
- Immutable raw data + mutable edit layer (non-destructive editing)
- Edit history stored as a list of (position, old_base, new_base) tuples for undo/redo
"""
from __future__ import annotations

import copy
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np


# Standard nucleotide colors used throughout the application
BASE_COLORS = {
    "A": "#00aa00",  # Green
    "T": "#cc0000",  # Red
    "G": "#000000",  # Black (visible in light mode; swapped in dark mode)
    "C": "#0000cc",  # Blue
}

# Dark-mode override for G (black → gold)
BASE_COLORS_DARK = {
    "A": "#33dd33",
    "T": "#ff4444",
    "G": "#dddd00",
    "C": "#4488ff",
}


@dataclass
class ChromatogramMetadata:
    """Instrument and run metadata extracted from file headers."""

    sample_name: str = ""
    file_path: str = ""
    file_format: str = ""  # "ab1" or "scf"
    instrument_model: str = ""
    run_date: str = ""
    lane: int = 0
    spacing: float = 0.0
    dye_set: str = ""
    comment: str = ""


@dataclass
class Chromatogram:
    """Central data model for a single chromatogram trace.

    Attributes:
        traces: dict mapping nucleotide letter ('A','T','G','C') to 1-D
            numpy array of fluorescence intensities (analyzed traces).
        raw_traces: same structure but for unprocessed DATA1-4 channels.
        basecalls: list of called nucleotide characters, one per peak.
        peak_locations: 1-D int array of sample-point indices where each
            base was called (indexes into the trace arrays).
        quality_scores: 1-D int array of Phred quality values per base.
        metadata: ChromatogramMetadata instance.
        trim_start: index into basecalls where good-quality region begins.
        trim_end: index into basecalls where good-quality region ends (exclusive).
    """

    traces: dict[str, np.ndarray] = field(default_factory=dict)
    raw_traces: dict[str, np.ndarray] = field(default_factory=dict)
    basecalls: list[str] = field(default_factory=list)
    peak_locations: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))
    quality_scores: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))
    metadata: ChromatogramMetadata = field(default_factory=ChromatogramMetadata)

    # Trimming boundaries (indices into basecalls)
    trim_start: int = 0
    trim_end: int = 0  # 0 means "use full length"

    # Edit tracking
    _edit_history: list[tuple[int, str, str]] = field(default_factory=list, repr=False)
    _redo_stack: list[tuple[int, str, str]] = field(default_factory=list, repr=False)
    _original_basecalls: list[str] = field(default_factory=list, repr=False)

    def __post_init__(self):
        if self.trim_end == 0 and len(self.basecalls) > 0:
            self.trim_end = len(self.basecalls)
        if not self._original_basecalls and self.basecalls:
            self._original_basecalls = list(self.basecalls)

    # ── Properties ──────────────────────────────────────────────────

    @property
    def num_bases(self) -> int:
        return len(self.basecalls)

    @property
    def trace_length(self) -> int:
        """Number of sample points in the longest trace channel."""
        if not self.traces:
            return 0
        return max(len(v) for v in self.traces.values())

    @property
    def trimmed_basecalls(self) -> list[str]:
        return self.basecalls[self.trim_start : self.trim_end]

    @property
    def trimmed_quality(self) -> np.ndarray:
        return self.quality_scores[self.trim_start : self.trim_end]

    @property
    def sequence(self) -> str:
        """Full called sequence as a string."""
        return "".join(self.basecalls)

    @property
    def trimmed_sequence(self) -> str:
        return "".join(self.trimmed_basecalls)

    @property
    def mean_quality(self) -> float:
        if len(self.quality_scores) == 0:
            return 0.0
        return float(np.mean(self.trimmed_quality))

    @property
    def is_edited(self) -> bool:
        return len(self._edit_history) > 0

    # ── Editing ─────────────────────────────────────────────────────

    def edit_base(self, position: int, new_base: str) -> None:
        """Change the base call at *position*, recording the edit for undo."""
        if position < 0 or position >= len(self.basecalls):
            raise IndexError(f"Position {position} out of range [0, {len(self.basecalls)})")
        new_base = new_base.upper()
        if new_base not in "ACGTNRYSWKMBDHV-":
            raise ValueError(f"Invalid base character: {new_base!r}")
        old_base = self.basecalls[position]
        if old_base == new_base:
            return
        self.basecalls[position] = new_base
        self._edit_history.append((position, old_base, new_base))
        self._redo_stack.clear()

    def undo(self) -> Optional[tuple[int, str, str]]:
        """Undo the last edit. Returns (position, old, new) or None."""
        if not self._edit_history:
            return None
        pos, old_base, new_base = self._edit_history.pop()
        self.basecalls[pos] = old_base
        self._redo_stack.append((pos, old_base, new_base))
        return (pos, old_base, new_base)

    def redo(self) -> Optional[tuple[int, str, str]]:
        """Redo the last undone edit."""
        if not self._redo_stack:
            return None
        pos, old_base, new_base = self._redo_stack.pop()
        self.basecalls[pos] = new_base
        self._edit_history.append((pos, old_base, new_base))
        return (pos, old_base, new_base)

    def reset_edits(self) -> None:
        """Revert all edits to original base calls."""
        if self._original_basecalls:
            self.basecalls = list(self._original_basecalls)
        self._edit_history.clear()
        self._redo_stack.clear()

    def get_edits(self) -> list[dict]:
        """Return a list of all edits as dicts for display / export."""
        edits = []
        for pos, old, new in self._edit_history:
            edits.append({"position": pos + 1, "original": old, "edited": new})
        return edits

    # ── Cloning ─────────────────────────────────────────────────────

    def clone(self) -> Chromatogram:
        """Deep copy of this chromatogram."""
        return copy.deepcopy(self)
