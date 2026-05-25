"""Quality score visualization widget.

Displays Phred quality scores as a bar chart synchronized with
the trace view. Color-coded: green (≥30), yellow (20-29), red (<20).
"""
from __future__ import annotations

from typing import Optional

import numpy as np

try:
    import pyqtgraph as pg
    from pyqtgraph import PlotWidget, InfiniteLine
    from PyQt6.QtWidgets import QWidget, QVBoxLayout
    HAS_GUI = True
except ImportError:
    HAS_GUI = False

from ..core.models import Chromatogram


if HAS_GUI:

    class QualityWidget(QWidget):
        """Bar chart of per-base Phred quality scores."""

        def __init__(self, parent=None, dark_mode: bool = True):
            super().__init__(parent)
            self._chrom: Optional[Chromatogram] = None
            self._dark_mode = dark_mode
            self._setup_ui()

        def _setup_ui(self):
            layout = QVBoxLayout(self)
            layout.setContentsMargins(0, 0, 0, 0)

            bg = "#1e1e2e" if self._dark_mode else "#fafafa"
            self._plot_widget = PlotWidget(background=bg)
            self._plot = self._plot_widget.getPlotItem()
            self._plot.setLabel("left", "Phred Score")
            self._plot.setLabel("bottom", "")
            self._plot.showGrid(x=False, y=True, alpha=0.15)
            self._plot.setMaximumHeight(120)

            # Q20 threshold line
            self._q20_line = InfiniteLine(
                pos=20, angle=0,
                pen=pg.mkPen("#ff8800", width=1, style=pg.QtCore.Qt.PenStyle.DashLine),
            )
            self._plot.addItem(self._q20_line)

            layout.addWidget(self._plot_widget)

        def set_chromatogram(self, chrom: Chromatogram):
            self._chrom = chrom
            self._plot.clear()
            self._plot.addItem(self._q20_line)

            if chrom.num_bases == 0:
                return

            locs = chrom.peak_locations.astype(float)
            quals = chrom.quality_scores.astype(float)

            # Color by quality level
            brushes = []
            for q in quals:
                if q >= 30:
                    brushes.append(pg.mkBrush("#2ecc71"))
                elif q >= 20:
                    brushes.append(pg.mkBrush("#f39c12"))
                else:
                    brushes.append(pg.mkBrush("#e74c3c"))

            bar = pg.BarGraphItem(
                x=locs, height=quals, width=6,
                brushes=brushes, pen=pg.mkPen(None),
            )
            self._plot.addItem(bar)

        def set_view_range(self, x_min: float, x_max: float):
            """Sync x-range with trace widget."""
            self._plot.setXRange(x_min, x_max, padding=0)

        def set_dark_mode(self, dark: bool):
            self._dark_mode = dark
            bg = "#1e1e2e" if dark else "#fafafa"
            self._plot_widget.setBackground(bg)
            if self._chrom:
                self.set_chromatogram(self._chrom)

else:
    class QualityWidget:
        def __init__(self, *args, **kwargs):
            raise ImportError("PyQt6 and pyqtgraph required for QualityWidget")
