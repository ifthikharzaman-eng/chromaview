"""Sidebar file explorer for chromatogram files.

Displays a tree view of a selected directory filtered to show
.ab1 and .scf files. Double-click opens a file.
"""
from __future__ import annotations

from pathlib import Path

try:
    from PyQt6.QtCore import Qt, pyqtSignal
    from PyQt6.QtWidgets import (
        QWidget, QVBoxLayout, QTreeView, QPushButton,
        QFileDialog, QLabel, QHBoxLayout,
    )
    from PyQt6.QtGui import QFileSystemModel, QIcon
    HAS_GUI = True
except ImportError:
    HAS_GUI = False


if HAS_GUI:

    class FileBrowser(QWidget):
        """Sidebar file browser filtered for chromatogram files.

        Signals:
            file_selected(str): Emitted when user double-clicks a .ab1/.scf file.
        """

        file_selected = pyqtSignal(str)

        SUPPORTED_EXTENSIONS = {".ab1", ".scf"}

        def __init__(self, parent=None):
            super().__init__(parent)
            self._current_dir = str(Path.home())
            self._setup_ui()

        def _setup_ui(self):
            layout = QVBoxLayout(self)
            layout.setContentsMargins(4, 4, 4, 4)

            # Header with folder picker
            header = QHBoxLayout()
            self._dir_label = QLabel("No folder selected")
            self._dir_label.setWordWrap(True)
            self._dir_label.setStyleSheet("font-size: 11px; color: #888;")
            header.addWidget(self._dir_label, stretch=1)

            browse_btn = QPushButton("📂")
            browse_btn.setToolTip("Choose folder")
            browse_btn.setFixedSize(30, 30)
            browse_btn.clicked.connect(self._choose_directory)
            header.addWidget(browse_btn)
            layout.addLayout(header)

            # File system model
            self._model = QFileSystemModel()
            self._model.setNameFilters(["*.ab1", "*.scf", "*.AB1", "*.SCF"])
            self._model.setNameFilterDisables(False)

            # Tree view
            self._tree = QTreeView()
            self._tree.setModel(self._model)
            self._tree.setRootIndex(self._model.setRootPath(self._current_dir))
            self._tree.setHeaderHidden(True)
            # Hide size, type, date columns
            for col in (1, 2, 3):
                self._tree.hideColumn(col)
            self._tree.doubleClicked.connect(self._on_double_click)
            self._tree.setAnimated(True)
            self._tree.setSortingEnabled(True)

            layout.addWidget(self._tree)

        def _choose_directory(self):
            directory = QFileDialog.getExistingDirectory(
                self, "Select Folder", self._current_dir,
            )
            if directory:
                self._current_dir = directory
                self._tree.setRootIndex(self._model.setRootPath(directory))
                self._dir_label.setText(directory)

        def set_directory(self, path: str):
            """Programmatically set the browsed directory."""
            self._current_dir = path
            self._tree.setRootIndex(self._model.setRootPath(path))
            self._dir_label.setText(path)

        def _on_double_click(self, index):
            filepath = self._model.filePath(index)
            if filepath and Path(filepath).suffix.lower() in self.SUPPORTED_EXTENSIONS:
                self.file_selected.emit(filepath)

else:
    class FileBrowser:
        def __init__(self, *args, **kwargs):
            raise ImportError("PyQt6 required for FileBrowser")
