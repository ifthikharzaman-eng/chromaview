"""BLAST dialog: BLASTWorker (QThread) + BLASTDialog + BLASTSettingsDialog."""
from __future__ import annotations

import io
from pathlib import Path
from typing import Optional

try:
    from PyQt6.QtCore import Qt, QSettings, QThread, pyqtSignal, QUrl
    from PyQt6.QtGui import QFont, QDesktopServices
    from PyQt6.QtWidgets import (
        QDialog, QVBoxLayout, QHBoxLayout, QGroupBox, QLabel,
        QPushButton, QComboBox, QLineEdit, QSpinBox, QDoubleSpinBox,
        QRadioButton, QButtonGroup, QFileDialog, QProgressBar,
        QTableWidget, QTableWidgetItem, QHeaderView, QSplitter,
        QTextEdit, QDialogButtonBox, QMessageBox, QWidget,
    )
    HAS_GUI = True
except ImportError:
    HAS_GUI = False

from ..analysis.blast import BLASTHit, DATABASES, MIN_QUERY_LENGTH, parse_blast_record
from ..core.models import Chromatogram

SETTINGS_EMAIL_KEY = "ncbi/email"
SETTINGS_TOOL_KEY = "ncbi/tool"
NCBI_NUCCORE_URL = "https://www.ncbi.nlm.nih.gov/nuccore/{accession}"


if HAS_GUI:

    class BLASTWorker(QThread):
        """Background thread that submits a blastn search and emits results."""

        progress = pyqtSignal(str)
        result_ready = pyqtSignal(list)   # list[BLASTHit]
        error = pyqtSignal(str)

        def __init__(
            self,
            sequence: str,
            database: str,
            email: str,
            tool: str = "ChromaView",
            hitlist_size: int = 20,
            expect: float = 10.0,
            parent=None,
        ):
            super().__init__(parent)
            self._sequence = sequence
            self._database = database
            self._email = email or None
            self._tool = tool or "ChromaView"
            self._hitlist_size = hitlist_size
            self._expect = expect
            self._cancelled = False

        def cancel(self) -> None:
            self._cancelled = True
            self.requestInterruption()

        def run(self) -> None:
            try:
                from Bio.Blast import NCBIWWW
                from Bio import Blast

                NCBIWWW.email = self._email
                NCBIWWW.tool = self._tool

                self.progress.emit("Submitting query to NCBI BLAST…")
                handle = NCBIWWW.qblast(
                    "blastn",
                    self._database,
                    self._sequence,
                    hitlist_size=self._hitlist_size,
                    expect=self._expect,
                    format_type="XML",
                )

                if self._cancelled:
                    return

                self.progress.emit("Parsing results…")
                # NCBIWWW.qblast returns StringIO (text mode); Bio.Blast.read
                # requires binary mode, so read the text and re-encode to bytes.
                raw = handle.read()
                if isinstance(raw, str):
                    raw = raw.encode("utf-8")
                record = Blast.read(io.BytesIO(raw))
                hits = parse_blast_record(record)
                self.result_ready.emit(hits)

            except Exception as exc:
                if not self._cancelled:
                    self.error.emit(str(exc))

    # ── Numeric table item (sorts by stored float, not display text) ──────

    class _NumItem(QTableWidgetItem):
        def __lt__(self, other: QTableWidgetItem) -> bool:
            try:
                return (
                    float(self.data(Qt.ItemDataRole.UserRole))
                    < float(other.data(Qt.ItemDataRole.UserRole))
                )
            except (TypeError, ValueError):
                return super().__lt__(other)

    # ── Settings dialog ───────────────────────────────────────────────────

    class BLASTSettingsDialog(QDialog):
        """Persistent NCBI email / tool-name settings."""

        def __init__(self, parent=None):
            super().__init__(parent)
            self.setWindowTitle("NCBI BLAST Settings")
            self.setMinimumWidth(380)

            s = QSettings()
            layout = QVBoxLayout(self)

            grp = QGroupBox("NCBI Etiquette Parameters")
            g = QVBoxLayout(grp)

            row1 = QHBoxLayout()
            row1.addWidget(QLabel("Contact email:"))
            self._email = QLineEdit(s.value(SETTINGS_EMAIL_KEY, ""))
            self._email.setPlaceholderText("you@example.com")
            row1.addWidget(self._email)
            g.addLayout(row1)

            row2 = QHBoxLayout()
            row2.addWidget(QLabel("Tool name:"))
            self._tool = QLineEdit(s.value(SETTINGS_TOOL_KEY, "ChromaView"))
            row2.addWidget(self._tool)
            g.addLayout(row2)

            note = QLabel(
                "NCBI requests that automated queries include a contact email "
                "so they can reach you if your script causes issues."
            )
            note.setWordWrap(True)
            note.setStyleSheet("color: #888; font-size: 11px;")
            g.addWidget(note)

            layout.addWidget(grp)

            btns = QDialogButtonBox(
                QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
            )
            btns.accepted.connect(self._save_and_accept)
            btns.rejected.connect(self.reject)
            layout.addWidget(btns)

        def _save_and_accept(self) -> None:
            s = QSettings()
            s.setValue(SETTINGS_EMAIL_KEY, self._email.text().strip())
            s.setValue(SETTINGS_TOOL_KEY, self._tool.text().strip() or "ChromaView")
            self.accept()

    # ── Main BLAST dialog ─────────────────────────────────────────────────

    class BLASTDialog(QDialog):
        """Run blastn against NCBI, show hits in a sortable table."""

        def __init__(self, parent=None, chrom: Optional[Chromatogram] = None):
            super().__init__(parent)
            self.setWindowTitle("NCBI BLAST Search")
            self.setMinimumSize(860, 680)
            self.resize(960, 740)

            self._chrom = chrom
            self._worker: Optional[BLASTWorker] = None
            self._hits: list[BLASTHit] = []

            layout = QVBoxLayout(self)
            layout.setSpacing(6)

            layout.addWidget(self._build_query_group())
            layout.addWidget(self._build_settings_group())
            layout.addLayout(self._build_run_row())
            layout.addWidget(self._build_results_splitter(), stretch=1)

            self._refresh_query_preview()

        # ── Query section ─────────────────────────────────────────

        def _build_query_group(self) -> QGroupBox:
            grp = QGroupBox("Query Sequence")
            g = QVBoxLayout(grp)

            self._radio_chrom = QRadioButton("Use current chromatogram (trimmed sequence)")
            self._radio_file = QRadioButton("Open a chromatogram file:")
            self._btn_grp = QButtonGroup(self)
            self._btn_grp.addButton(self._radio_chrom, 0)
            self._btn_grp.addButton(self._radio_file, 1)

            if self._chrom:
                self._radio_chrom.setChecked(True)
            else:
                self._radio_chrom.setEnabled(False)
                self._radio_file.setChecked(True)

            g.addWidget(self._radio_chrom)

            file_row = QHBoxLayout()
            file_row.addWidget(self._radio_file)
            self._file_path_label = QLabel("(no file selected)")
            self._file_path_label.setStyleSheet("color: #888;")
            file_row.addWidget(self._file_path_label, stretch=1)
            self._browse_btn = QPushButton("Browse…")
            self._browse_btn.clicked.connect(self._browse_file)
            file_row.addWidget(self._browse_btn)
            g.addLayout(file_row)

            self._seq_preview = QLabel("No sequence loaded.")
            self._seq_preview.setStyleSheet("font-family: monospace; font-size: 11px; color: #aaa;")
            self._seq_preview.setWordWrap(True)
            g.addWidget(self._seq_preview)

            self._radio_chrom.toggled.connect(self._refresh_query_preview)
            self._radio_file.toggled.connect(self._refresh_query_preview)

            return grp

        def _build_settings_group(self) -> QGroupBox:
            grp = QGroupBox("Search Settings")
            g = QHBoxLayout(grp)

            s = QSettings()

            g.addWidget(QLabel("Database:"))
            self._db_combo = QComboBox()
            for db_id, db_label in DATABASES:
                self._db_combo.addItem(db_label, db_id)
            self._db_combo.setCurrentIndex(0)
            g.addWidget(self._db_combo)

            g.addSpacing(12)
            g.addWidget(QLabel("Max hits:"))
            self._hits_spin = QSpinBox()
            self._hits_spin.setRange(1, 100)
            self._hits_spin.setValue(20)
            g.addWidget(self._hits_spin)

            g.addSpacing(12)
            g.addWidget(QLabel("E-value:"))
            self._expect_spin = QDoubleSpinBox()
            self._expect_spin.setRange(1e-200, 1000.0)
            self._expect_spin.setValue(10.0)
            self._expect_spin.setDecimals(1)
            g.addWidget(self._expect_spin)

            g.addSpacing(12)
            g.addWidget(QLabel("NCBI email:"))
            self._email_edit = QLineEdit(s.value(SETTINGS_EMAIL_KEY, ""))
            self._email_edit.setPlaceholderText("you@example.com (recommended)")
            self._email_edit.setMinimumWidth(180)
            g.addWidget(self._email_edit)

            g.addStretch()
            return grp

        def _build_run_row(self) -> QHBoxLayout:
            row = QHBoxLayout()

            self._run_btn = QPushButton("Run BLAST")
            self._run_btn.setDefault(True)
            self._run_btn.clicked.connect(self._run_blast)
            row.addWidget(self._run_btn)

            self._cancel_btn = QPushButton("Cancel Search")
            self._cancel_btn.setEnabled(False)
            self._cancel_btn.clicked.connect(self._cancel_search)
            row.addWidget(self._cancel_btn)

            row.addSpacing(12)
            self._status_label = QLabel("Ready")
            row.addWidget(self._status_label, stretch=1)

            self._progress = QProgressBar()
            self._progress.setRange(0, 0)  # indeterminate
            self._progress.setVisible(False)
            self._progress.setMaximumWidth(180)
            row.addWidget(self._progress)

            close_btn = QPushButton("Close")
            close_btn.clicked.connect(self.reject)
            row.addWidget(close_btn)

            return row

        def _build_results_splitter(self) -> QSplitter:
            splitter = QSplitter(Qt.Orientation.Vertical)

            # ── Top: hits table ──────────────────────────────────
            top = QWidget()
            top_layout = QVBoxLayout(top)
            top_layout.setContentsMargins(0, 0, 0, 0)

            self._hits_label = QLabel("Results: (no search run yet)")
            top_layout.addWidget(self._hits_label)

            self._table = QTableWidget(0, 7)
            self._table.setHorizontalHeaderLabels([
                "Accession", "Description", "Organism",
                "Query Cov%", "Identity%", "E-value", "Bit Score",
            ])
            self._table.setSortingEnabled(True)
            self._table.setSelectionBehavior(
                QTableWidget.SelectionBehavior.SelectRows
            )
            self._table.setEditTriggers(
                QTableWidget.EditTrigger.NoEditTriggers
            )
            hdr = self._table.horizontalHeader()
            hdr.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
            hdr.setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
            hdr.setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
            for col in (3, 4, 5, 6):
                hdr.setSectionResizeMode(col, QHeaderView.ResizeMode.ResizeToContents)
            self._table.currentCellChanged.connect(self._on_hit_selected)
            top_layout.addWidget(self._table)
            splitter.addWidget(top)

            # ── Bottom: alignment viewer ──────────────────────────
            bottom = QWidget()
            bot_layout = QVBoxLayout(bottom)
            bot_layout.setContentsMargins(0, 0, 0, 0)

            bot_layout.addWidget(QLabel("Alignment:"))
            self._align_view = QTextEdit()
            self._align_view.setReadOnly(True)
            mono = QFont("Courier New", 10)
            mono.setStyleHint(QFont.StyleHint.Monospace)
            self._align_view.setFont(mono)
            bot_layout.addWidget(self._align_view)

            ncbi_row = QHBoxLayout()
            self._ncbi_btn = QPushButton("Open in NCBI")
            self._ncbi_btn.setEnabled(False)
            self._ncbi_btn.clicked.connect(self._open_in_ncbi)
            ncbi_row.addWidget(self._ncbi_btn)
            ncbi_row.addStretch()
            bot_layout.addLayout(ncbi_row)

            splitter.addWidget(bottom)
            splitter.setSizes([420, 220])
            return splitter

        # ── Query preview ─────────────────────────────────────────

        def _refresh_query_preview(self) -> None:
            seq = self._get_query_sequence()
            if seq:
                n = len(seq)
                preview = seq[:70] + ("…" if n > 70 else "")
                warn = f"  ⚠ Too short (<{MIN_QUERY_LENGTH} bp)" if n < MIN_QUERY_LENGTH else ""
                self._seq_preview.setText(f"{preview}  [{n} bp]{warn}")
            else:
                self._seq_preview.setText("No sequence loaded.")

        def _get_query_sequence(self) -> str:
            if self._radio_chrom.isChecked() and self._chrom:
                return self._chrom.trimmed_sequence
            if self._radio_file.isChecked() and hasattr(self, "_file_chrom") and self._file_chrom:
                return self._file_chrom.trimmed_sequence
            return ""

        def _browse_file(self) -> None:
            filepath, _ = QFileDialog.getOpenFileName(
                self, "Open Chromatogram for BLAST",
                "", "Chromatogram files (*.ab1 *.scf);;All (*)",
            )
            if not filepath:
                return
            try:
                p = Path(filepath)
                from ..core.ab1_parser import parse_ab1
                from ..core.scf_parser import parse_scf
                from ..analysis.quality import trim_low_quality

                if p.suffix.lower() == ".ab1":
                    chrom = parse_ab1(p)
                elif p.suffix.lower() == ".scf":
                    chrom = parse_scf(p)
                else:
                    QMessageBox.warning(self, "Unsupported file", f"Cannot read {p.suffix}")
                    return

                trim_low_quality(chrom, threshold=20)
                self._file_chrom = chrom
                self._file_path_label.setText(p.name)
                self._file_path_label.setStyleSheet("")
                self._radio_file.setChecked(True)
                self._refresh_query_preview()
            except Exception as exc:
                QMessageBox.critical(self, "Error loading file", str(exc))

        # ── BLAST run / cancel ────────────────────────────────────

        def _run_blast(self) -> None:
            seq = self._get_query_sequence()
            if not seq:
                QMessageBox.warning(self, "No sequence", "Please load a chromatogram first.")
                return

            if len(seq) < MIN_QUERY_LENGTH:
                ans = QMessageBox.question(
                    self, "Short sequence",
                    f"Sequence is only {len(seq)} bp (< {MIN_QUERY_LENGTH} bp recommended). "
                    "BLAST results may be unreliable. Continue?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                )
                if ans != QMessageBox.StandardButton.Yes:
                    return

            email = self._email_edit.text().strip()
            QSettings().setValue(SETTINGS_EMAIL_KEY, email)

            db = self._db_combo.currentData()
            self._worker = BLASTWorker(
                sequence=seq,
                database=db,
                email=email,
                hitlist_size=self._hits_spin.value(),
                expect=self._expect_spin.value(),
                parent=self,
            )
            self._worker.progress.connect(self._on_progress)
            self._worker.result_ready.connect(self._on_results)
            self._worker.error.connect(self._on_error)
            self._worker.finished.connect(self._on_finished)

            self._run_btn.setEnabled(False)
            self._cancel_btn.setEnabled(True)
            self._progress.setVisible(True)
            self._table.setRowCount(0)
            self._align_view.clear()
            self._ncbi_btn.setEnabled(False)
            self._hits = []
            self._worker.start()

        def _cancel_search(self) -> None:
            if self._worker and self._worker.isRunning():
                self._worker.cancel()
                self._worker.wait(3000)
            self._status_label.setText("Cancelled.")
            self._on_finished()

        # ── Worker signal handlers ────────────────────────────────

        def _on_progress(self, msg: str) -> None:
            self._status_label.setText(msg)

        def _on_results(self, hits: list) -> None:
            self._hits = hits
            self._populate_table(hits)
            n = len(hits)
            self._hits_label.setText(f"Results: {n} hit{'s' if n != 1 else ''}")
            if n == 0:
                self._status_label.setText("Search complete — no hits found.")
                self._align_view.setPlainText(
                    "No significant alignments found.\n\n"
                    "Tips:\n"
                    "• Try a different database (e.g. 'nt')\n"
                    "• Increase the E-value cutoff\n"
                    "• Check that your sequence is correct"
                )
            else:
                self._status_label.setText(f"Done — {n} hit{'s' if n != 1 else ''} found.")

        def _on_error(self, msg: str) -> None:
            self._status_label.setText(f"Error: {msg[:80]}")
            err_text = msg
            if "urlopen" in msg.lower() or "connection" in msg.lower() or "network" in msg.lower():
                err_text = (
                    "Network error — could not reach NCBI BLAST.\n\n"
                    "Please check your internet connection and try again.\n\n"
                    f"Details: {msg}"
                )
            elif "timeout" in msg.lower():
                err_text = (
                    "NCBI BLAST timed out.\n\n"
                    "BLAST searches can take 30–120 seconds. "
                    "Please try again later or during off-peak hours.\n\n"
                    f"Details: {msg}"
                )
            QMessageBox.critical(self, "BLAST Error", err_text)

        def _on_finished(self) -> None:
            self._run_btn.setEnabled(True)
            self._cancel_btn.setEnabled(False)
            self._progress.setVisible(False)

        # ── Results table ─────────────────────────────────────────

        def _populate_table(self, hits: list[BLASTHit]) -> None:
            self._table.setSortingEnabled(False)
            self._table.setRowCount(0)

            for hit in hits:
                row = self._table.rowCount()
                self._table.insertRow(row)

                self._table.setItem(row, 0, QTableWidgetItem(hit.accession))
                # Truncate description for display
                desc = hit.description if len(hit.description) <= 60 else hit.description[:57] + "…"
                self._table.setItem(row, 1, QTableWidgetItem(desc))
                self._table.setItem(row, 2, QTableWidgetItem(hit.organism))

                cov_item = _NumItem(f"{hit.query_coverage:.1%}")
                cov_item.setData(Qt.ItemDataRole.UserRole, hit.query_coverage)
                self._table.setItem(row, 3, cov_item)

                id_item = _NumItem(f"{hit.pct_identity:.1%}")
                id_item.setData(Qt.ItemDataRole.UserRole, hit.pct_identity)
                self._table.setItem(row, 4, id_item)

                eval_str = f"{hit.evalue:.2e}" if hit.evalue < 0.01 else f"{hit.evalue:.4f}"
                eval_item = _NumItem(eval_str)
                eval_item.setData(Qt.ItemDataRole.UserRole, hit.evalue)
                self._table.setItem(row, 5, eval_item)

                bs_item = _NumItem(f"{hit.bit_score:.1f}")
                bs_item.setData(Qt.ItemDataRole.UserRole, hit.bit_score)
                self._table.setItem(row, 6, bs_item)

            self._table.setSortingEnabled(True)

        def _on_hit_selected(self, row: int, _col: int, _old_row: int, _old_col: int) -> None:
            if row < 0 or row >= len(self._hits):
                self._align_view.clear()
                self._ncbi_btn.setEnabled(False)
                return
            # Find the hit corresponding to this table row
            # (table may be sorted, so we need to match by accession)
            acc_item = self._table.item(row, 0)
            if acc_item is None:
                return
            acc = acc_item.text()
            hit = next((h for h in self._hits if h.accession == acc), None)
            if hit is None:
                return
            self._align_view.setPlainText(hit.alignment_text)
            self._ncbi_btn.setEnabled(True)
            self._ncbi_btn.setProperty("accession", acc)

        def _open_in_ncbi(self) -> None:
            acc = self._ncbi_btn.property("accession") or ""
            if acc:
                url = QUrl(NCBI_NUCCORE_URL.format(accession=acc))
                QDesktopServices.openUrl(url)

        def reject(self) -> None:
            if self._worker and self._worker.isRunning():
                self._worker.cancel()
                self._worker.wait(3000)
            super().reject()

    # ── No-GUI stubs ──────────────────────────────────────────────────────

else:
    class BLASTWorker:
        def __init__(self, *a, **kw):
            raise ImportError("PyQt6 required")

    class BLASTDialog:
        def __init__(self, *a, **kw):
            raise ImportError("PyQt6 required")

    class BLASTSettingsDialog:
        def __init__(self, *a, **kw):
            raise ImportError("PyQt6 required")
