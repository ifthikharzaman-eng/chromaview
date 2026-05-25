"""Main application window.

Assembles all GUI components into a single window with:
- Menu bar (File, Edit, Analysis, View, Help)
- Toolbar (open, zoom, theme toggle, search)
- Central area: chromatogram trace + sequence bar
- Dockable panels: file browser (left), quality chart (bottom)
- Status bar with sample info
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

try:
    from PyQt6.QtCore import Qt, QSettings
    from PyQt6.QtGui import QAction, QKeySequence, QIcon, QFont
    from PyQt6.QtWidgets import (
        QMainWindow, QWidget, QVBoxLayout, QDockWidget,
        QFileDialog, QMessageBox, QStatusBar, QToolBar,
        QLineEdit, QLabel, QApplication, QHBoxLayout,
    )
    HAS_GUI = True
except ImportError:
    HAS_GUI = False

from ..core.models import Chromatogram
from ..core.ab1_parser import parse_ab1
from ..core.scf_parser import parse_scf
from ..core.sequence_ops import reverse_complement, align_to_reference, search_sequence
from ..analysis.quality import trim_low_quality, quality_summary
from ..analysis.peaks import average_snr
from ..analysis.mutation import detect_ambiguous_bases, call_mutations
from ..export.fasta import export_fasta
from ..export.csv_export import export_csv
from ..export.image import export_image


if HAS_GUI:

    from .trace_widget import TraceWidget
    from .sequence_bar import SequenceBar
    from .quality_widget import QualityWidget
    from .file_browser import FileBrowser
    from .dialogs import TrimDialog, CompareDialog, BatchDialog
    from .blast_dialog import BLASTDialog, BLASTSettingsDialog
    from .consensus_dialog import ConsensusDialog
    from .theme import DARK_THEME, LIGHT_THEME, get_stylesheet, ThemeColors

    class MainWindow(QMainWindow):
        """ChromaView main application window."""

        APP_NAME = "ChromaView"
        APP_VERSION = "0.1.0"

        def __init__(self):
            super().__init__()
            self._chrom: Optional[Chromatogram] = None
            self._dark_mode = True
            self._theme = DARK_THEME

            self.setWindowTitle(self.APP_NAME)
            self.setMinimumSize(1100, 700)
            self.resize(1400, 850)

            self._create_central()
            self._create_docks()
            self._create_menus()
            self._create_toolbar()
            self._create_statusbar()
            self._apply_theme()
            self._connect_signals()

        # ── Central Widget ───────────────────────────────────────

        def _create_central(self):
            central = QWidget()
            layout = QVBoxLayout(central)
            layout.setContentsMargins(0, 0, 0, 0)
            layout.setSpacing(0)

            self._trace_widget = TraceWidget(dark_mode=self._dark_mode)
            layout.addWidget(self._trace_widget, stretch=1)

            self._sequence_bar = SequenceBar(dark_mode=self._dark_mode)
            layout.addWidget(self._sequence_bar)

            self.setCentralWidget(central)

        # ── Dock Widgets ─────────────────────────────────────────

        def _create_docks(self):
            # File browser dock (left)
            self._file_browser = FileBrowser()
            file_dock = QDockWidget("File Browser", self)
            file_dock.setWidget(self._file_browser)
            file_dock.setFeatures(
                QDockWidget.DockWidgetFeature.DockWidgetMovable |
                QDockWidget.DockWidgetFeature.DockWidgetClosable
            )
            self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, file_dock)
            self._file_dock = file_dock

            # Quality chart dock (bottom)
            self._quality_widget = QualityWidget(dark_mode=self._dark_mode)
            quality_dock = QDockWidget("Quality Scores", self)
            quality_dock.setWidget(self._quality_widget)
            quality_dock.setFeatures(
                QDockWidget.DockWidgetFeature.DockWidgetMovable |
                QDockWidget.DockWidgetFeature.DockWidgetClosable
            )
            self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, quality_dock)
            self._quality_dock = quality_dock

            # Info panel dock (right)
            self._info_panel = QWidget()
            info_layout = QVBoxLayout(self._info_panel)
            self._info_text = QLabel("No file loaded")
            self._info_text.setWordWrap(True)
            self._info_text.setAlignment(Qt.AlignmentFlag.AlignTop)
            self._info_text.setStyleSheet("padding: 8px; font-size: 12px;")
            info_layout.addWidget(self._info_text)
            info_layout.addStretch()

            info_dock = QDockWidget("Info", self)
            info_dock.setWidget(self._info_panel)
            info_dock.setFeatures(
                QDockWidget.DockWidgetFeature.DockWidgetMovable |
                QDockWidget.DockWidgetFeature.DockWidgetClosable
            )
            self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, info_dock)
            self._info_dock = info_dock

        # ── Menus ────────────────────────────────────────────────

        def _create_menus(self):
            menubar = self.menuBar()

            # ── File ──
            file_menu = menubar.addMenu("&File")

            open_act = QAction("&Open...", self)
            open_act.setShortcut(QKeySequence.StandardKey.Open)
            open_act.triggered.connect(self._open_file)
            file_menu.addAction(open_act)

            file_menu.addSeparator()

            export_fasta_act = QAction("Export &FASTA...", self)
            export_fasta_act.setShortcut(QKeySequence("Ctrl+Shift+F"))
            export_fasta_act.triggered.connect(self._export_fasta)
            file_menu.addAction(export_fasta_act)

            export_csv_act = QAction("Export &CSV...", self)
            export_csv_act.triggered.connect(self._export_csv)
            file_menu.addAction(export_csv_act)

            export_img_act = QAction("Export &Image...", self)
            export_img_act.setShortcut(QKeySequence("Ctrl+Shift+I"))
            export_img_act.triggered.connect(self._export_image)
            file_menu.addAction(export_img_act)

            file_menu.addSeparator()

            quit_act = QAction("&Quit", self)
            quit_act.setShortcut(QKeySequence.StandardKey.Quit)
            quit_act.triggered.connect(self.close)
            file_menu.addAction(quit_act)

            # ── Edit ──
            edit_menu = menubar.addMenu("&Edit")

            undo_act = QAction("&Undo", self)
            undo_act.setShortcut(QKeySequence.StandardKey.Undo)
            undo_act.triggered.connect(self._undo)
            edit_menu.addAction(undo_act)

            redo_act = QAction("&Redo", self)
            redo_act.setShortcut(QKeySequence.StandardKey.Redo)
            redo_act.triggered.connect(self._redo)
            edit_menu.addAction(redo_act)

            edit_menu.addSeparator()

            reset_act = QAction("Reset All Edits", self)
            reset_act.triggered.connect(self._reset_edits)
            edit_menu.addAction(reset_act)

            # ── Analysis ──
            analysis_menu = menubar.addMenu("&Analysis")

            trim_act = QAction("&Trim Low-Quality Ends...", self)
            trim_act.setShortcut(QKeySequence("Ctrl+T"))
            trim_act.triggered.connect(self._trim_dialog)
            analysis_menu.addAction(trim_act)

            revcomp_act = QAction("&Reverse Complement", self)
            revcomp_act.setShortcut(QKeySequence("Ctrl+R"))
            revcomp_act.triggered.connect(self._reverse_complement)
            analysis_menu.addAction(revcomp_act)

            compare_act = QAction("&Compare to Reference...", self)
            compare_act.setShortcut(QKeySequence("Ctrl+Shift+C"))
            compare_act.triggered.connect(self._compare_dialog)
            analysis_menu.addAction(compare_act)

            ambig_act = QAction("Detect &Ambiguous Bases", self)
            ambig_act.triggered.connect(self._detect_ambiguous)
            analysis_menu.addAction(ambig_act)

            analysis_menu.addSeparator()

            consensus_act = QAction("Build &Consensus (F+R)...", self)
            consensus_act.setShortcut(QKeySequence("Ctrl+Shift+N"))
            consensus_act.triggered.connect(self._build_consensus)
            analysis_menu.addAction(consensus_act)

            analysis_menu.addSeparator()

            self._blast_act = QAction("Run &BLAST...", self)
            self._blast_act.setShortcut(QKeySequence("Ctrl+B"))
            self._blast_act.setEnabled(False)
            self._blast_act.triggered.connect(self._run_blast)
            analysis_menu.addAction(self._blast_act)

            blast_file_act = QAction("BLAST from &File...", self)
            blast_file_act.triggered.connect(self._blast_from_file)
            analysis_menu.addAction(blast_file_act)

            analysis_menu.addSeparator()

            batch_act = QAction("&Batch Processing...", self)
            batch_act.triggered.connect(self._batch_dialog)
            analysis_menu.addAction(batch_act)

            # ── View ──
            view_menu = menubar.addMenu("&View")

            theme_act = QAction("Toggle &Dark/Light Mode", self)
            theme_act.setShortcut(QKeySequence("Ctrl+D"))
            theme_act.triggered.connect(self._toggle_theme)
            view_menu.addAction(theme_act)

            view_menu.addSeparator()
            view_menu.addAction(self._file_dock.toggleViewAction())
            view_menu.addAction(self._quality_dock.toggleViewAction())
            view_menu.addAction(self._info_dock.toggleViewAction())

            zoom_in_act = QAction("Zoom &In", self)
            zoom_in_act.setShortcut(QKeySequence.StandardKey.ZoomIn)
            zoom_in_act.triggered.connect(self._zoom_in)
            view_menu.addAction(zoom_in_act)

            zoom_out_act = QAction("Zoom &Out", self)
            zoom_out_act.setShortcut(QKeySequence.StandardKey.ZoomOut)
            zoom_out_act.triggered.connect(self._zoom_out)
            view_menu.addAction(zoom_out_act)

            fit_act = QAction("&Fit to Window", self)
            fit_act.setShortcut(QKeySequence("Ctrl+0"))
            fit_act.triggered.connect(self._fit_view)
            view_menu.addAction(fit_act)

            # ── Help ──
            help_menu = menubar.addMenu("&Help")
            about_act = QAction("&About ChromaView", self)
            about_act.triggered.connect(self._show_about)
            help_menu.addAction(about_act)

            ncbi_settings_act = QAction("NCBI BLAST &Settings...", self)
            ncbi_settings_act.triggered.connect(self._ncbi_settings)
            help_menu.addAction(ncbi_settings_act)

        # ── Toolbar ──────────────────────────────────────────────

        def _create_toolbar(self):
            toolbar = QToolBar("Main Toolbar")
            toolbar.setMovable(False)
            self.addToolBar(toolbar)

            open_btn = QAction("📂 Open", self)
            open_btn.triggered.connect(self._open_file)
            toolbar.addAction(open_btn)

            toolbar.addSeparator()

            zoom_in = QAction("🔍+ Zoom In", self)
            zoom_in.triggered.connect(self._zoom_in)
            toolbar.addAction(zoom_in)

            zoom_out = QAction("🔍− Zoom Out", self)
            zoom_out.triggered.connect(self._zoom_out)
            toolbar.addAction(zoom_out)

            fit_btn = QAction("⬜ Fit", self)
            fit_btn.triggered.connect(self._fit_view)
            toolbar.addAction(fit_btn)

            toolbar.addSeparator()

            theme_btn = QAction("🌓 Theme", self)
            theme_btn.triggered.connect(self._toggle_theme)
            toolbar.addAction(theme_btn)

            toolbar.addSeparator()

            self._blast_toolbar_act = QAction("🧬 BLAST", self)
            self._blast_toolbar_act.setEnabled(False)
            self._blast_toolbar_act.setToolTip("Run NCBI BLAST on the loaded chromatogram (Ctrl+B)")
            self._blast_toolbar_act.triggered.connect(self._run_blast)
            toolbar.addAction(self._blast_toolbar_act)

            toolbar.addSeparator()

            # Search field
            self._search_field = QLineEdit()
            self._search_field.setPlaceholderText("Search position or motif...")
            self._search_field.setMaximumWidth(200)
            self._search_field.returnPressed.connect(self._search)
            toolbar.addWidget(self._search_field)

            search_btn = QAction("🔎", self)
            search_btn.triggered.connect(self._search)
            toolbar.addAction(search_btn)

        # ── Status Bar ───────────────────────────────────────────

        def _create_statusbar(self):
            self._status = QStatusBar()
            self.setStatusBar(self._status)
            self._status_label = QLabel("Ready")
            self._status_label.setObjectName("status_label")
            self._status.addWidget(self._status_label)

        # ── Signals ──────────────────────────────────────────────

        def _connect_signals(self):
            self._file_browser.file_selected.connect(self._load_file)
            self._trace_widget.base_selected.connect(self._on_base_selected)
            self._trace_widget.base_edited.connect(self._on_base_edited)
            self._trace_widget.view_range_changed.connect(self._on_view_range_changed)
            self._sequence_bar.base_clicked.connect(self._on_base_clicked)

        # ── Theme ────────────────────────────────────────────────

        def _apply_theme(self):
            self._theme = DARK_THEME if self._dark_mode else LIGHT_THEME
            self.setStyleSheet(get_stylesheet(self._theme))

        def _toggle_theme(self):
            self._dark_mode = not self._dark_mode
            self._apply_theme()
            self._trace_widget.set_dark_mode(self._dark_mode)
            self._sequence_bar.set_dark_mode(self._dark_mode)
            self._quality_widget.set_dark_mode(self._dark_mode)

        # ── File Operations ──────────────────────────────────────

        def _open_file(self):
            filepath, _ = QFileDialog.getOpenFileName(
                self, "Open Chromatogram",
                "", "Chromatogram files (*.ab1 *.scf);;AB1 files (*.ab1);;SCF files (*.scf);;All (*)",
            )
            if filepath:
                self._load_file(filepath)

        def _load_file(self, filepath: str):
            try:
                p = Path(filepath)
                if p.suffix.lower() == ".ab1":
                    chrom = parse_ab1(p)
                elif p.suffix.lower() == ".scf":
                    chrom = parse_scf(p)
                else:
                    QMessageBox.warning(self, "Error", f"Unsupported file format: {p.suffix}")
                    return

                # Auto-trim
                trim_low_quality(chrom, threshold=20)

                self._chrom = chrom
                self._display_chromatogram()
                self._status_label.setText(
                    f"{p.name}  |  {chrom.num_bases} bases  |  "
                    f"Mean Q: {chrom.mean_quality:.1f}  |  "
                    f"Trimmed: {chrom.trim_start}-{chrom.trim_end}"
                )

                # Set file browser to same directory
                self._file_browser.set_directory(str(p.parent))

            except Exception as e:
                QMessageBox.critical(self, "Error Loading File", str(e))

        def _display_chromatogram(self):
            if not self._chrom:
                return
            self._trace_widget.set_chromatogram(self._chrom)
            self._sequence_bar.set_chromatogram(self._chrom)
            self._quality_widget.set_chromatogram(self._chrom)
            self._update_info()
            self._blast_act.setEnabled(True)
            self._blast_toolbar_act.setEnabled(True)

        def _update_info(self):
            if not self._chrom:
                return
            qs = quality_summary(self._chrom)
            snr = average_snr(self._chrom)
            ambig = detect_ambiguous_bases(self._chrom)

            info = f"""<b>{self._chrom.metadata.sample_name or 'Unknown'}</b><br>
<b>File:</b> {Path(self._chrom.metadata.file_path).name}<br>
<b>Format:</b> {self._chrom.metadata.file_format.upper()}<br>
<b>Instrument:</b> {self._chrom.metadata.instrument_model or 'N/A'}<br>
<hr>
<b>Total bases:</b> {qs.total_bases}<br>
<b>Trimmed region:</b> {self._chrom.trim_start + 1}–{self._chrom.trim_end}<br>
<b>Usable bases:</b> {qs.trimmed_bases}<br>
<hr>
<b>Mean quality:</b> {qs.mean_quality:.1f}<br>
<b>Median quality:</b> {qs.median_quality:.1f}<br>
<b>Q≥20:</b> {qs.q20_count} ({qs.q20_fraction:.1%})<br>
<b>Q≥30:</b> {qs.q30_count} ({qs.q30_fraction:.1%})<br>
<b>Median SNR:</b> {snr:.1f}<br>
<hr>
<b>Ambiguous bases:</b> {len(ambig)}<br>
<b>Edits:</b> {len(self._chrom.get_edits())}
"""
            self._info_text.setText(info)

        # ── Signal Handlers ──────────────────────────────────────

        def _on_base_selected(self, index: int):
            self._sequence_bar.set_selected_base(index)
            if self._chrom and index < len(self._chrom.quality_scores):
                q = self._chrom.quality_scores[index]
                base = self._chrom.basecalls[index]
                self._status_label.setText(
                    f"Position {index + 1}: {base}  |  Quality: {q}"
                )

        def _on_base_edited(self, index: int, new_base: str):
            self._sequence_bar.set_chromatogram(self._chrom)
            self._update_info()

        def _on_view_range_changed(self, x_min: float, x_max: float):
            self._sequence_bar.set_view_range(x_min, x_max)
            self._quality_widget.set_view_range(x_min, x_max)

        def _on_base_clicked(self, index: int):
            self._trace_widget.scroll_to_base(index)
            self._on_base_selected(index)

        # ── Edit Actions ─────────────────────────────────────────

        def _undo(self):
            if self._chrom:
                result = self._chrom.undo()
                if result:
                    self._display_chromatogram()

        def _redo(self):
            if self._chrom:
                result = self._chrom.redo()
                if result:
                    self._display_chromatogram()

        def _reset_edits(self):
            if self._chrom and self._chrom.is_edited:
                self._chrom.reset_edits()
                self._display_chromatogram()

        # ── Analysis Actions ─────────────────────────────────────

        def _trim_dialog(self):
            if not self._chrom:
                return
            dlg = TrimDialog(self, current_threshold=20)
            if dlg.exec():
                trim_low_quality(self._chrom, threshold=dlg.threshold, window=dlg.window)
                self._display_chromatogram()
                self._status_label.setText(
                    f"Trimmed: {self._chrom.trim_start + 1}–{self._chrom.trim_end}"
                )

        def _reverse_complement(self):
            if not self._chrom:
                return
            seq = self._chrom.sequence
            rc = reverse_complement(seq)
            QMessageBox.information(
                self, "Reverse Complement",
                f"Original length: {len(seq)}\n\n"
                f"Reverse complement:\n{rc[:100]}{'...' if len(rc) > 100 else ''}",
            )

        def _compare_dialog(self):
            if not self._chrom:
                return
            dlg = CompareDialog(self)
            if dlg.exec():
                ref_seq = dlg.get_sequence()
                if not ref_seq:
                    return
                result = align_to_reference(self._chrom.trimmed_sequence, ref_seq)
                mutations = call_mutations(
                    self._chrom.trimmed_sequence, ref_seq,
                    quality_scores=self._chrom.trimmed_quality,
                    alignment_result=result,
                )

                # Map mismatch positions back to full-sequence indices
                mismatch_full = [p + self._chrom.trim_start for p in result.mismatch_positions]
                self._trace_widget.highlight_bases(mismatch_full)
                self._sequence_bar.set_highlights(mismatch_full)

                msg = (
                    f"Identity: {result.identity:.1%}\n"
                    f"Score: {result.score}\n"
                    f"Mismatches: {len(result.mismatch_positions)}\n"
                    f"Mutations found: {len(mutations)}\n\n"
                )
                for m in mutations[:20]:
                    msg += f"  Pos {m.query_position + 1}: {m.ref_base}→{m.query_base} ({m.mutation_type}, Q={m.quality})\n"
                if len(mutations) > 20:
                    msg += f"  ... and {len(mutations) - 20} more\n"

                QMessageBox.information(self, "Comparison Result", msg)

        def _detect_ambiguous(self):
            if not self._chrom:
                return
            ambigs = detect_ambiguous_bases(self._chrom)
            if not ambigs:
                QMessageBox.information(self, "Ambiguous Bases", "No ambiguous bases detected.")
                return
            positions = [a.position for a in ambigs]
            self._trace_widget.highlight_bases(positions, color="#ff8800")
            self._sequence_bar.set_highlights(positions)

            msg = f"Found {len(ambigs)} ambiguous positions:\n\n"
            for a in ambigs[:30]:
                msg += (
                    f"  Pos {a.position + 1}: {a.called_base}→{a.suggested_iupac} "
                    f"(secondary {a.secondary_base}, ratio {a.ratio:.2f})\n"
                )
            QMessageBox.information(self, "Ambiguous Bases", msg)

        def _run_blast(self):
            dlg = BLASTDialog(self, chrom=self._chrom)
            dlg.exec()

        def _blast_from_file(self):
            dlg = BLASTDialog(self, chrom=None)
            dlg._radio_file.setChecked(True)
            dlg._browse_btn.click()
            dlg.exec()

        def _ncbi_settings(self):
            BLASTSettingsDialog(self).exec()

        def _build_consensus(self):
            dlg = ConsensusDialog(self, chrom=self._chrom)
            dlg.exec()

        def _batch_dialog(self):
            dlg = BatchDialog(self)
            if dlg.exec():
                files = dlg.files
                if not files:
                    return

                output_dir = QFileDialog.getExistingDirectory(
                    self, "Select Output Directory"
                )
                if not output_dir:
                    return

                dlg.progress.setVisible(True)
                dlg.progress.setMaximum(len(files))

                for i, filepath in enumerate(files):
                    try:
                        p = Path(filepath)
                        if p.suffix.lower() == ".ab1":
                            chrom = parse_ab1(p)
                        else:
                            chrom = parse_scf(p)

                        trim_low_quality(chrom, threshold=dlg.threshold_spin.value())
                        qs = quality_summary(chrom)

                        if dlg.export_fasta_cb.isChecked():
                            export_fasta(chrom, Path(output_dir) / f"{p.stem}.fasta")
                        if dlg.export_csv_cb.isChecked():
                            export_csv(chrom, Path(output_dir) / f"{p.stem}.csv")
                        if dlg.export_image_cb.isChecked():
                            export_image(chrom, Path(output_dir) / f"{p.stem}.png")

                        dlg.add_result(
                            p.name, qs.trimmed_bases, qs.mean_quality,
                            qs.q20_fraction, "✓ OK",
                        )
                    except Exception as e:
                        dlg.add_result(p.name, 0, 0, 0, f"✗ {e}")

                    dlg.progress.setValue(i + 1)
                    QApplication.processEvents()

                dlg.progress.setVisible(False)
                self._status_label.setText(f"Batch complete: {len(files)} files processed")

        # ── Search ───────────────────────────────────────────────

        def _search(self):
            if not self._chrom:
                return
            query = self._search_field.text().strip()
            if not query:
                return
            hits = search_sequence(self._chrom.sequence, query)
            if hits:
                positions = [h.start for h in hits]
                self._trace_widget.highlight_bases(positions, color="#00aaff")
                self._sequence_bar.set_highlights(positions)
                self._trace_widget.scroll_to_base(positions[0])
                self._status_label.setText(f"Found {len(hits)} match(es) for '{query}'")
            else:
                self._status_label.setText(f"No matches for '{query}'")

        # ── View Actions ─────────────────────────────────────────

        def _zoom_in(self):
            vb = self._trace_widget._plot.getViewBox()
            vb.scaleBy((0.5, 1))

        def _zoom_out(self):
            vb = self._trace_widget._plot.getViewBox()
            vb.scaleBy((2, 1))

        def _fit_view(self):
            self._trace_widget._plot.autoRange()

        # ── Export Actions ───────────────────────────────────────

        def _export_fasta(self):
            if not self._chrom:
                return
            filepath, _ = QFileDialog.getSaveFileName(
                self, "Export FASTA", "", "FASTA (*.fasta *.fa);;All (*)",
            )
            if filepath:
                export_fasta(self._chrom, filepath)
                self._status_label.setText(f"Exported: {filepath}")

        def _export_csv(self):
            if not self._chrom:
                return
            filepath, _ = QFileDialog.getSaveFileName(
                self, "Export CSV", "", "CSV (*.csv);;All (*)",
            )
            if filepath:
                export_csv(self._chrom, filepath, include_traces=True)
                self._status_label.setText(f"Exported: {filepath}")

        def _export_image(self):
            if not self._chrom:
                return
            filepath, _ = QFileDialog.getSaveFileName(
                self, "Export Image", "",
                "PNG (*.png);;SVG (*.svg);;All (*)",
            )
            if filepath:
                try:
                    export_image(
                        self._chrom, filepath,
                        title=self._chrom.metadata.sample_name,
                        show_quality=True,
                    )
                    self._status_label.setText(f"Exported: {filepath}")
                except ImportError as e:
                    QMessageBox.warning(self, "Export Error", str(e))

        # ── About ────────────────────────────────────────────────

        def _show_about(self):
            QMessageBox.about(
                self, f"About {self.APP_NAME}",
                f"<b>{self.APP_NAME} v{self.APP_VERSION}</b><br><br>"
                "Open-source Sanger sequencing chromatogram analyzer.<br><br>"
                "Built with Python, PyQt6, PyQtGraph, and Biopython.<br><br>"
                "© 2026 ChromaView Contributors<br>"
                "MIT License",
            )

        # ── Drag and Drop ────────────────────────────────────────

        def dragEnterEvent(self, event):
            if event.mimeData().hasUrls():
                event.acceptProposedAction()

        def dropEvent(self, event):
            for url in event.mimeData().urls():
                filepath = url.toLocalFile()
                if Path(filepath).suffix.lower() in (".ab1", ".scf"):
                    self._load_file(filepath)
                    break

else:
    class MainWindow:
        def __init__(self, *args, **kwargs):
            raise ImportError("PyQt6 required for MainWindow")
