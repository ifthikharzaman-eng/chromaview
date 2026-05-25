"""Consensus assembly dialog for paired forward + reverse Sanger reads."""
from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

try:
    from PyQt6.QtCore import Qt, QSettings
    from PyQt6.QtGui import QFont
    from PyQt6.QtWidgets import (
        QDialog, QVBoxLayout, QHBoxLayout, QGroupBox, QLabel,
        QPushButton, QLineEdit, QFileDialog, QSplitter,
        QTextEdit, QMessageBox, QWidget,
    )
    HAS_GUI = True
except ImportError:
    HAS_GUI = False

from ..analysis.consensus import (
    ConsensusResult,
    build_consensus,
    find_overlap,
    reverse_complement_with_quality,
)
from ..core.models import Chromatogram, ChromatogramMetadata

# ── filename pair-detection patterns ─────────────────────────────────────────
# Each tuple: (forward_token, reverse_token) — applied case-insensitively to
# the filename stem (with surrounding non-alpha boundary).

_PAIR_TOKENS = [
    ("_F_",       "_R_"),
    ("-F_",       "-R_"),
    ("_fwd_",     "_rev_"),
    ("_forward_", "_reverse_"),
]
_PAIR_STEMS = [   # tokens that may appear at the very end of the stem
    ("_F",        "_R"),
    ("-F",        "-R"),
    ("_fwd",      "_rev"),
    ("_forward",  "_reverse"),
]


def _suggest_pair(path: Path) -> Optional[Path]:
    """Try to find the paired read filename by swapping F↔R tokens."""
    name = path.name
    stem = path.stem
    ext  = path.suffix
    directory = path.parent

    for fwd_tok, rev_tok in _PAIR_TOKENS:
        pattern = re.compile(re.escape(fwd_tok), re.IGNORECASE)
        candidate = pattern.sub(
            lambda m: _match_case(rev_tok, m.group()), name, count=1
        )
        if candidate != name:
            p = directory / candidate
            if p.exists():
                return p

    for fwd_stem, rev_stem in _PAIR_STEMS:
        pattern = re.compile(re.escape(fwd_stem) + r"$", re.IGNORECASE)
        new_stem = pattern.sub(
            lambda m: _match_case(rev_stem, m.group()), stem, count=1
        )
        if new_stem != stem:
            p = directory / (new_stem + ext)
            if p.exists():
                return p

    return None


def _match_case(template: str, original: str) -> str:
    """Return template in the same case style as original."""
    if original.isupper():
        return template.upper()
    if original.islower():
        return template.lower()
    if original[1:].islower():   # Title-like: _Fwd
        return template[0] + template[1:].lower()
    return template


if HAS_GUI:
    import numpy as np

    def _load_chrom(path: Path) -> Chromatogram:
        """Parse AB1 or SCF and auto-trim."""
        from ..core.ab1_parser import parse_ab1
        from ..core.scf_parser import parse_scf
        from ..analysis.quality import trim_low_quality

        if path.suffix.lower() == ".ab1":
            chrom = parse_ab1(path)
        elif path.suffix.lower() == ".scf":
            chrom = parse_scf(path)
        else:
            raise ValueError(f"Unsupported format: {path.suffix}")
        trim_low_quality(chrom, threshold=20)
        return chrom

    # ── Alignment formatter ───────────────────────────────────────────────────

    def _format_alignment(result: ConsensusResult, cols: int = 60) -> str:
        """Return a human-readable alignment report."""
        ov = result.overlap
        lines: list[str] = []

        if not ov.has_overlap:
            lines.append("⚠ No overlap detected between the two reads.")
            lines.append("")
            lines.append("Possible causes:")
            lines.append("  • Reads are not from the same amplicon")
            lines.append("  • One read is in the wrong orientation (try Swap ↕)")
            lines.append("  • Very low quality / heavily trimmed reads")
            return "\n".join(lines)

        # ── Stats ─────────────────────────────────────────────────────────
        lines.append(
            f"Overlap: {ov.overlap_length} bp  |  "
            f"Identity: {ov.identity:.1%}  |  "
            f"Disagreements: {result.disagreements}  |  "
            f"Ambiguities: {result.ambiguities}"
        )
        lines.append(
            f"Consensus: {result.total_length} bp  "
            f"(fwd-only: {result.fwd_only_length}  "
            f"overlap: {result.overlap_length}  "
            f"rev-only: {result.rev_only_length})"
        )
        lines.append("")

        # ── Warnings ──────────────────────────────────────────────────────
        for w in result.warnings:
            lines.append(f"⚠ {w}")
        if result.warnings:
            lines.append("")

        # ── Wrapped pairwise alignment (overlap region only) ───────────────
        lines.append("─── Overlap alignment ───")
        fwd_aln = ov.fwd_aligned
        rc_aln  = ov.rc_aligned
        for start in range(0, len(fwd_aln), cols):
            fa = fwd_aln[start : start + cols]
            ra = rc_aln [start : start + cols]
            mid = "".join(
                "|" if f != "-" and r != "-" and f.upper() == r.upper() else " "
                for f, r in zip(fa, ra)
            )
            q_offset = ov.fwd_start + sum(c != "-" for c in fwd_aln[:start])
            lines.append(f"fwd {q_offset + 1:>5}  {fa}")
            lines.append(f"       {'':>5}  {mid}")
            q_offset2 = ov.rc_start + sum(c != "-" for c in rc_aln[:start])
            lines.append(f"rc  {q_offset2 + 1:>5}  {ra}")
            lines.append("")

        # ── Consensus preview ──────────────────────────────────────────────
        lines.append("─── Consensus (first 120 bp) ───")
        seq = result.consensus_seq
        for start in range(0, min(120, len(seq)), cols):
            lines.append(seq[start : start + cols])
        if len(seq) > 120:
            lines.append(f"... ({len(seq) - 120} more bases)")

        return "\n".join(lines)

    # ── Dialog ────────────────────────────────────────────────────────────────

    class ConsensusDialog(QDialog):
        """Build and inspect a consensus from a forward + reverse read pair."""

        def __init__(self, parent=None, chrom: Optional[Chromatogram] = None):
            super().__init__(parent)
            self.setWindowTitle("Build Consensus — Forward + Reverse Reads")
            self.setMinimumSize(820, 640)
            self.resize(920, 720)

            self._fwd_chrom: Optional[Chromatogram] = chrom
            self._fwd_path: Optional[Path] = (
                Path(chrom.metadata.file_path) if chrom and chrom.metadata.file_path else None
            )
            self._rev_chrom: Optional[Chromatogram] = None
            self._rev_path: Optional[Path] = None
            self._result: Optional[ConsensusResult] = None

            layout = QVBoxLayout(self)
            layout.setSpacing(6)
            layout.addWidget(self._build_pair_group())
            layout.addLayout(self._build_action_row())
            layout.addWidget(self._build_results_widget(), stretch=1)

            # Pre-populate and auto-pair
            if self._fwd_path:
                self._fwd_path_edit.setText(str(self._fwd_path))
                self._try_auto_pair()

        # ── Pair selection ────────────────────────────────────────────────

        def _build_pair_group(self) -> QGroupBox:
            grp = QGroupBox("Read Pair")
            layout = QVBoxLayout(grp)

            # Forward row
            fwd_row = QHBoxLayout()
            fwd_row.addWidget(QLabel("Forward read:"))
            self._fwd_path_edit = QLineEdit()
            self._fwd_path_edit.setPlaceholderText("Path to forward .ab1/.scf file")
            self._fwd_path_edit.setReadOnly(True)
            fwd_row.addWidget(self._fwd_path_edit, stretch=1)
            fwd_browse = QPushButton("Browse…")
            fwd_browse.clicked.connect(self._browse_fwd)
            fwd_row.addWidget(fwd_browse)
            layout.addLayout(fwd_row)

            # Reverse row
            rev_row = QHBoxLayout()
            rev_row.addWidget(QLabel("Reverse read:"))
            self._rev_path_edit = QLineEdit()
            self._rev_path_edit.setPlaceholderText("Path to reverse .ab1/.scf file")
            self._rev_path_edit.setReadOnly(True)
            rev_row.addWidget(self._rev_path_edit, stretch=1)
            rev_browse = QPushButton("Browse…")
            rev_browse.clicked.connect(self._browse_rev)
            rev_row.addWidget(rev_browse)
            layout.addLayout(rev_row)

            # Helper buttons
            btn_row = QHBoxLayout()
            auto_btn = QPushButton("Auto-pair")
            auto_btn.setToolTip(
                "Detect the reverse read filename by swapping F↔R tokens "
                "(_F/_R, _fwd/_rev, _forward/_reverse)"
            )
            auto_btn.clicked.connect(self._try_auto_pair)
            btn_row.addWidget(auto_btn)

            swap_btn = QPushButton("Swap ↕")
            swap_btn.setToolTip("Swap forward and reverse assignments")
            swap_btn.clicked.connect(self._swap_reads)
            btn_row.addWidget(swap_btn)
            btn_row.addStretch()

            note = QLabel(
                "The reverse read will be reverse-complemented automatically before alignment."
            )
            note.setStyleSheet("color: #888; font-size: 11px;")
            btn_row.addWidget(note)
            layout.addLayout(btn_row)

            return grp

        def _build_action_row(self) -> QHBoxLayout:
            row = QHBoxLayout()
            self._build_btn = QPushButton("Build Consensus")
            self._build_btn.setDefault(True)
            self._build_btn.clicked.connect(self._build_consensus)
            row.addWidget(self._build_btn)

            self._status_label = QLabel("")
            row.addWidget(self._status_label, stretch=1)

            close_btn = QPushButton("Close")
            close_btn.clicked.connect(self.reject)
            row.addWidget(close_btn)
            return row

        def _build_results_widget(self) -> QWidget:
            container = QWidget()
            layout = QVBoxLayout(container)
            layout.setContentsMargins(0, 0, 0, 0)

            self._result_view = QTextEdit()
            self._result_view.setReadOnly(True)
            mono = QFont("Courier New", 10)
            mono.setStyleHint(QFont.StyleHint.Monospace)
            self._result_view.setFont(mono)
            self._result_view.setPlaceholderText(
                "Load a forward and reverse read pair, then click Build Consensus."
            )
            layout.addWidget(self._result_view)

            # Export / BLAST buttons
            btn_row = QHBoxLayout()
            self._blast_btn = QPushButton("Send to BLAST")
            self._blast_btn.setEnabled(False)
            self._blast_btn.clicked.connect(self._send_to_blast)
            btn_row.addWidget(self._blast_btn)

            self._export_btn = QPushButton("Export FASTA…")
            self._export_btn.setEnabled(False)
            self._export_btn.clicked.connect(self._export_fasta)
            btn_row.addWidget(self._export_btn)

            self._copy_btn = QPushButton("Copy Sequence")
            self._copy_btn.setEnabled(False)
            self._copy_btn.clicked.connect(self._copy_sequence)
            btn_row.addWidget(self._copy_btn)

            btn_row.addStretch()
            layout.addLayout(btn_row)
            return container

        # ── File browsing ─────────────────────────────────────────────────

        def _browse_fwd(self) -> None:
            path = self._browse_chrom("Select Forward Chromatogram")
            if path:
                self._fwd_path = path
                self._fwd_chrom = None  # will be re-loaded on Build
                self._fwd_path_edit.setText(str(path))

        def _browse_rev(self) -> None:
            path = self._browse_chrom("Select Reverse Chromatogram")
            if path:
                self._rev_path = path
                self._rev_chrom = None
                self._rev_path_edit.setText(str(path))

        def _browse_chrom(self, title: str) -> Optional[Path]:
            filepath, _ = QFileDialog.getOpenFileName(
                self, title, "",
                "Chromatogram files (*.ab1 *.scf);;AB1 (*.ab1);;SCF (*.scf);;All (*)",
            )
            return Path(filepath) if filepath else None

        def _try_auto_pair(self) -> None:
            if not self._fwd_path:
                return
            suggested = _suggest_pair(self._fwd_path)
            if suggested:
                self._rev_path = suggested
                self._rev_chrom = None
                self._rev_path_edit.setText(str(suggested))
                self._status_label.setText(f"Auto-paired: {suggested.name}")
            else:
                self._status_label.setText("No matching reverse read found automatically.")

        def _swap_reads(self) -> None:
            self._fwd_path, self._rev_path = self._rev_path, self._fwd_path
            self._fwd_chrom, self._rev_chrom = self._rev_chrom, self._fwd_chrom
            self._fwd_path_edit.setText(str(self._fwd_path) if self._fwd_path else "")
            self._rev_path_edit.setText(str(self._rev_path) if self._rev_path else "")

        # ── Consensus building ────────────────────────────────────────────

        def _build_consensus(self) -> None:
            if not self._fwd_path or not self._rev_path:
                QMessageBox.warning(
                    self, "Missing reads",
                    "Please select both a forward and a reverse read file."
                )
                return

            try:
                self._status_label.setText("Loading files…")
                fwd = _load_chrom(self._fwd_path)
                rev = _load_chrom(self._rev_path)
            except Exception as exc:
                QMessageBox.critical(self, "Error loading files", str(exc))
                self._status_label.setText("Error loading files.")
                return

            fwd_seq = fwd.trimmed_sequence
            fwd_q   = np.array(fwd.trimmed_quality, dtype=int)
            rev_seq = rev.trimmed_sequence
            rev_q   = np.array(rev.trimmed_quality, dtype=int)

            if len(fwd_seq) < 20:
                QMessageBox.warning(
                    self, "Forward read too short",
                    f"Forward trimmed sequence is only {len(fwd_seq)} bp. "
                    "Check the file and trim settings."
                )
                return
            if len(rev_seq) < 20:
                QMessageBox.warning(
                    self, "Reverse read too short",
                    f"Reverse trimmed sequence is only {len(rev_seq)} bp."
                )
                return

            self._status_label.setText("Reverse-complementing reverse read…")
            rc_seq, rc_q = reverse_complement_with_quality(rev_seq, rev_q)

            self._status_label.setText("Aligning (Smith-Waterman)…")
            try:
                overlap = find_overlap(fwd_seq, rc_seq)
            except Exception as exc:
                QMessageBox.critical(self, "Alignment error", str(exc))
                self._status_label.setText("Alignment failed.")
                return

            if not overlap.has_overlap:
                self._result = None
                identity_line = (
                    f"\n\nBest alignment found: {overlap.overlap_length} bp "
                    f"at {overlap.identity:.1%} identity "
                    f"(minimum required: 90%).\n"
                    "This indicates a spurious match — the reads are likely "
                    "from different samples or organisms."
                    if overlap.overlap_length > 0 else ""
                )
                self._result_view.setPlainText(
                    "No valid overlap found between these two reads."
                    + identity_line
                    + "\n\nPossible causes:\n"
                    "  - Reads are not from the same sample (different well numbers?)\n"
                    "  - Reads are not from the same amplicon\n"
                    "  - Reads are in the same orientation "
                    "(try Swap to try the other direction)\n"
                    "  - One or both reads are very short / low quality"
                )
                self._blast_btn.setEnabled(False)
                self._export_btn.setEnabled(False)
                self._copy_btn.setEnabled(False)
                self._status_label.setText(
                    f"No valid overlap — best identity {overlap.identity:.1%} "
                    f"(need >= 90%)."
                    if overlap.overlap_length > 0
                    else "No overlap detected."
                )
                return

            self._status_label.setText("Building consensus…")
            self._result = build_consensus(fwd_seq, fwd_q, rc_seq, rc_q, overlap)

            text = _format_alignment(self._result)
            self._result_view.setPlainText(text)
            self._blast_btn.setEnabled(True)
            self._export_btn.setEnabled(True)
            self._copy_btn.setEnabled(True)

            ov = self._result.overlap
            self._status_label.setText(
                f"Done — {self._result.total_length} bp consensus  |  "
                f"overlap {ov.overlap_length} bp @ {ov.identity:.1%} identity"
            )

        # ── Export / BLAST ────────────────────────────────────────────────

        def _send_to_blast(self) -> None:
            if not self._result:
                return
            consensus_chrom = Chromatogram(
                basecalls=list(self._result.consensus_seq),
                quality_scores=np.array(self._result.consensus_quality, dtype=int),
                metadata=ChromatogramMetadata(sample_name="Consensus"),
            )
            from .blast_dialog import BLASTDialog
            dlg = BLASTDialog(self.parent(), chrom=consensus_chrom)
            dlg.exec()

        def _export_fasta(self) -> None:
            if not self._result:
                return
            fwd_stem = self._fwd_path.stem if self._fwd_path else "fwd"
            rev_stem = self._rev_path.stem if self._rev_path else "rev"
            default  = f"consensus_{fwd_stem}+{rev_stem}.fasta"
            filepath, _ = QFileDialog.getSaveFileName(
                self, "Export Consensus FASTA", default,
                "FASTA (*.fasta *.fa);;All (*)",
            )
            if not filepath:
                return

            ov = self._result.overlap
            header = (
                f"consensus {fwd_stem}+{rev_stem} "
                f"len:{self._result.total_length} "
                f"overlap:{ov.overlap_length} "
                f"identity:{ov.identity:.3f} "
                f"disagreements:{self._result.disagreements} "
                f"ambiguities:{self._result.ambiguities}"
            )
            seq = self._result.consensus_seq
            lines = [f">{header}"]
            for i in range(0, len(seq), 80):
                lines.append(seq[i : i + 80])
            Path(filepath).write_text("\n".join(lines) + "\n", encoding="utf-8")
            self._status_label.setText(f"Exported: {Path(filepath).name}")

        def _copy_sequence(self) -> None:
            if not self._result:
                return
            from PyQt6.QtWidgets import QApplication
            QApplication.clipboard().setText(self._result.consensus_seq)
            self._status_label.setText("Consensus sequence copied to clipboard.")

    # ── No-GUI stub ───────────────────────────────────────────────────────────

else:
    class ConsensusDialog:
        def __init__(self, *a, **kw):
            raise ImportError("PyQt6 required for ConsensusDialog")
