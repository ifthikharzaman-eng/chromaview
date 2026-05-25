"""Dialog windows for trim settings, reference comparison, and batch processing."""
from __future__ import annotations

from pathlib import Path
from typing import Optional

try:
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import (
        QDialog, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
        QSpinBox, QGroupBox, QTextEdit, QFileDialog, QProgressBar,
        QTableWidget, QTableWidgetItem, QHeaderView, QDialogButtonBox,
        QCheckBox,
    )
    HAS_GUI = True
except ImportError:
    HAS_GUI = False


if HAS_GUI:

    class TrimDialog(QDialog):
        """Dialog for configuring quality-based trimming parameters."""

        def __init__(self, parent=None, current_threshold: int = 20):
            super().__init__(parent)
            self.setWindowTitle("Trim Low-Quality Ends")
            self.setMinimumWidth(350)

            layout = QVBoxLayout(self)

            # Threshold setting
            group = QGroupBox("Trimming Parameters")
            g_layout = QVBoxLayout(group)

            row = QHBoxLayout()
            row.addWidget(QLabel("Quality threshold (Phred):"))
            self.threshold_spin = QSpinBox()
            self.threshold_spin.setRange(0, 60)
            self.threshold_spin.setValue(current_threshold)
            row.addWidget(self.threshold_spin)
            g_layout.addLayout(row)

            row2 = QHBoxLayout()
            row2.addWidget(QLabel("Minimum window size:"))
            self.window_spin = QSpinBox()
            self.window_spin.setRange(1, 50)
            self.window_spin.setValue(10)
            row2.addWidget(self.window_spin)
            g_layout.addLayout(row2)

            layout.addWidget(group)

            info = QLabel(
                "The modified Mott algorithm trims regions where the running "
                "quality score falls below the threshold. A Phred score of 20 "
                "corresponds to 99% base-call accuracy."
            )
            info.setWordWrap(True)
            info.setStyleSheet("color: #888; font-size: 11px;")
            layout.addWidget(info)

            # Buttons
            buttons = QDialogButtonBox(
                QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
            )
            buttons.accepted.connect(self.accept)
            buttons.rejected.connect(self.reject)
            layout.addWidget(buttons)

        @property
        def threshold(self) -> int:
            return self.threshold_spin.value()

        @property
        def window(self) -> int:
            return self.window_spin.value()


    class CompareDialog(QDialog):
        """Dialog for entering a reference sequence to compare against."""

        def __init__(self, parent=None):
            super().__init__(parent)
            self.setWindowTitle("Compare to Reference Sequence")
            self.setMinimumSize(500, 400)

            layout = QVBoxLayout(self)

            layout.addWidget(QLabel("Paste reference sequence (FASTA or raw):"))

            self.text_edit = QTextEdit()
            self.text_edit.setPlaceholderText(
                ">reference_sequence\nATGCGATCGATCGA..."
            )
            self.text_edit.setFont(
                self.text_edit.font()
            )
            layout.addWidget(self.text_edit)

            # Or load from file
            file_row = QHBoxLayout()
            file_row.addWidget(QLabel("Or load from file:"))
            load_btn = QPushButton("Browse...")
            load_btn.clicked.connect(self._load_file)
            file_row.addWidget(load_btn)
            layout.addLayout(file_row)

            buttons = QDialogButtonBox(
                QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
            )
            buttons.accepted.connect(self.accept)
            buttons.rejected.connect(self.reject)
            layout.addWidget(buttons)

        def _load_file(self):
            filepath, _ = QFileDialog.getOpenFileName(
                self, "Open Reference",
                "", "FASTA files (*.fasta *.fa *.fna);;Text files (*.txt);;All (*)",
            )
            if filepath:
                text = Path(filepath).read_text()
                self.text_edit.setPlainText(text)

        def get_sequence(self) -> str:
            """Extract the raw sequence from the text edit content."""
            text = self.text_edit.toPlainText().strip()
            lines = text.split("\n")
            seq_lines = [l.strip() for l in lines if not l.startswith(">")]
            return "".join(seq_lines).replace(" ", "").upper()


    class BatchDialog(QDialog):
        """Dialog for batch processing multiple chromatogram files."""

        def __init__(self, parent=None):
            super().__init__(parent)
            self.setWindowTitle("Batch Processing")
            self.setMinimumSize(700, 500)
            self._files: list[str] = []

            layout = QVBoxLayout(self)

            # File selection
            file_row = QHBoxLayout()
            add_btn = QPushButton("Add Files...")
            add_btn.clicked.connect(self._add_files)
            file_row.addWidget(add_btn)
            add_dir_btn = QPushButton("Add Folder...")
            add_dir_btn.clicked.connect(self._add_folder)
            file_row.addWidget(add_dir_btn)
            self._file_count_label = QLabel("0 files selected")
            file_row.addWidget(self._file_count_label)
            file_row.addStretch()
            layout.addLayout(file_row)

            # Options
            opts = QGroupBox("Export Options")
            opts_layout = QVBoxLayout(opts)
            self.export_fasta_cb = QCheckBox("Export trimmed FASTA")
            self.export_fasta_cb.setChecked(True)
            opts_layout.addWidget(self.export_fasta_cb)
            self.export_csv_cb = QCheckBox("Export quality CSV")
            self.export_csv_cb.setChecked(True)
            opts_layout.addWidget(self.export_csv_cb)
            self.export_image_cb = QCheckBox("Export chromatogram images")
            opts_layout.addWidget(self.export_image_cb)

            trim_row = QHBoxLayout()
            trim_row.addWidget(QLabel("Quality threshold:"))
            self.threshold_spin = QSpinBox()
            self.threshold_spin.setRange(0, 60)
            self.threshold_spin.setValue(20)
            trim_row.addWidget(self.threshold_spin)
            trim_row.addStretch()
            opts_layout.addLayout(trim_row)

            layout.addWidget(opts)

            # Results table
            self.results_table = QTableWidget()
            self.results_table.setColumnCount(5)
            self.results_table.setHorizontalHeaderLabels([
                "File", "Bases", "Mean Q", "Q20%", "Status",
            ])
            self.results_table.horizontalHeader().setSectionResizeMode(
                0, QHeaderView.ResizeMode.Stretch
            )
            layout.addWidget(self.results_table)

            # Progress
            self.progress = QProgressBar()
            self.progress.setVisible(False)
            layout.addWidget(self.progress)

            # Buttons
            btn_row = QHBoxLayout()
            self.run_btn = QPushButton("Run Batch")
            self.run_btn.clicked.connect(self.accept)
            btn_row.addWidget(self.run_btn)
            cancel_btn = QPushButton("Close")
            cancel_btn.clicked.connect(self.reject)
            btn_row.addWidget(cancel_btn)
            layout.addLayout(btn_row)

        def _add_files(self):
            files, _ = QFileDialog.getOpenFileNames(
                self, "Select Chromatogram Files",
                "", "Chromatogram files (*.ab1 *.scf);;All (*)",
            )
            self._files.extend(files)
            self._file_count_label.setText(f"{len(self._files)} files selected")

        def _add_folder(self):
            folder = QFileDialog.getExistingDirectory(self, "Select Folder")
            if folder:
                p = Path(folder)
                for ext in ("*.ab1", "*.scf", "*.AB1", "*.SCF"):
                    self._files.extend(str(f) for f in p.glob(ext))
                self._file_count_label.setText(f"{len(self._files)} files selected")

        @property
        def files(self) -> list[str]:
            return list(set(self._files))

        def add_result(self, filename: str, bases: int, mean_q: float,
                       q20_pct: float, status: str):
            """Add a row to the results table."""
            row = self.results_table.rowCount()
            self.results_table.insertRow(row)
            self.results_table.setItem(row, 0, QTableWidgetItem(filename))
            self.results_table.setItem(row, 1, QTableWidgetItem(str(bases)))
            self.results_table.setItem(row, 2, QTableWidgetItem(f"{mean_q:.1f}"))
            self.results_table.setItem(row, 3, QTableWidgetItem(f"{q20_pct:.1%}"))
            self.results_table.setItem(row, 4, QTableWidgetItem(status))


else:
    class TrimDialog:
        def __init__(self, *args, **kwargs):
            raise ImportError("PyQt6 required")

    class CompareDialog:
        def __init__(self, *args, **kwargs):
            raise ImportError("PyQt6 required")

    class BatchDialog:
        def __init__(self, *args, **kwargs):
            raise ImportError("PyQt6 required")
