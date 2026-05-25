"""Interactive chromatogram trace visualization widget.

Displays the four fluorescence channels (A/T/G/C) as colored curves,
with called base letters positioned at their peak locations. Supports
smooth zoom, pan, and scroll. Peak positions are marked with vertical
tick marks connecting the base letter to the trace.

Design decisions:
- PyQtGraph for fast GPU-accelerated 2D plotting with 10,000+ data points
- Custom TextItem overlay for base letters (cheaper than individual labels)
- Trim regions shown as semi-transparent shaded areas
- Editing: double-click a base to open an inline editor
"""
from __future__ import annotations

from typing import Optional

import numpy as np

try:
    import pyqtgraph as pg
    from pyqtgraph import PlotWidget, InfiniteLine, LinearRegionItem, TextItem
    from PyQt6.QtCore import Qt, pyqtSignal, QPointF
    from PyQt6.QtGui import QFont, QColor, QPen
    from PyQt6.QtWidgets import QVBoxLayout, QWidget, QInputDialog
    HAS_GUI = True
except ImportError:
    HAS_GUI = False

from ..core.models import Chromatogram, BASE_COLORS, BASE_COLORS_DARK


if HAS_GUI:

    class TraceWidget(QWidget):
        """Main chromatogram trace display widget.

        Signals:
            base_selected(int): Emitted when user clicks near a base.
            base_edited(int, str): Emitted after a base is edited.
            view_range_changed(float, float): Emitted on zoom/pan.
        """

        base_selected = pyqtSignal(int)
        base_edited = pyqtSignal(int, str)
        view_range_changed = pyqtSignal(float, float)

        def __init__(self, parent=None, dark_mode: bool = True):
            super().__init__(parent)
            self._chrom: Optional[Chromatogram] = None
            self._dark_mode = dark_mode
            self._base_items: list[TextItem] = []
            self._trace_curves: dict[str, pg.PlotDataItem] = {}
            self._trim_left: Optional[LinearRegionItem] = None
            self._trim_right: Optional[LinearRegionItem] = None
            self._cursor_line: Optional[InfiniteLine] = None
            self._highlight_items: list = []

            self._setup_ui()

        def _setup_ui(self):
            layout = QVBoxLayout(self)
            layout.setContentsMargins(0, 0, 0, 0)

            # Configure PyQtGraph
            bg = "#1e1e2e" if self._dark_mode else "#fafafa"
            self._plot_widget = PlotWidget(background=bg)
            self._plot = self._plot_widget.getPlotItem()
            self._plot.showGrid(x=False, y=False)
            self._plot.setLabel("left", "Fluorescence")
            self._plot.setLabel("bottom", "Sample Point")
            self._plot.getViewBox().setMouseMode(pg.ViewBox.RectMode)

            # Enable scroll-to-zoom on x-axis only
            self._plot.getViewBox().setMouseEnabled(x=True, y=True)
            self._plot.getViewBox().sigRangeChanged.connect(self._on_range_changed)

            # Cursor line
            self._cursor_line = InfiniteLine(
                angle=90, movable=False,
                pen=pg.mkPen("#ffffff44" if self._dark_mode else "#00000022", width=1)
            )
            self._plot.addItem(self._cursor_line)

            layout.addWidget(self._plot_widget)

            # Click handling
            self._plot_widget.scene().sigMouseClicked.connect(self._on_click)

        def set_chromatogram(self, chrom: Chromatogram):
            """Load and display a chromatogram."""
            self._chrom = chrom
            self._clear_display()
            self._draw_traces()
            self._draw_bases()
            self._draw_trim_regions()
            self._auto_range()

        def _clear_display(self):
            """Remove all plot items."""
            for curve in self._trace_curves.values():
                self._plot.removeItem(curve)
            self._trace_curves.clear()

            for item in self._base_items:
                self._plot.removeItem(item)
            self._base_items.clear()

            for item in self._highlight_items:
                self._plot.removeItem(item)
            self._highlight_items.clear()

            if self._trim_left:
                self._plot.removeItem(self._trim_left)
                self._trim_left = None
            if self._trim_right:
                self._plot.removeItem(self._trim_right)
                self._trim_right = None

        def _draw_traces(self):
            """Plot the four fluorescence channels."""
            if not self._chrom:
                return

            colors = BASE_COLORS_DARK if self._dark_mode else BASE_COLORS

            for base in ["A", "T", "G", "C"]:
                trace = self._chrom.traces.get(base)
                if trace is None or len(trace) == 0:
                    continue

                color = colors.get(base, "#888888")
                pen = pg.mkPen(color, width=1.2)
                x = np.arange(len(trace))
                curve = self._plot.plot(x, trace, pen=pen, name=base)
                self._trace_curves[base] = curve

        def _draw_bases(self):
            """Place base call letters at peak positions."""
            if not self._chrom:
                return

            colors = BASE_COLORS_DARK if self._dark_mode else BASE_COLORS
            font = QFont("Consolas", 9, QFont.Weight.Bold)

            # Find baseline for text placement (below the traces)
            y_offset = -50  # Place below x-axis

            for i, (base, loc) in enumerate(
                zip(self._chrom.basecalls, self._chrom.peak_locations)
            ):
                color = colors.get(base.upper(), "#888888")
                item = TextItem(text=base, color=color, anchor=(0.5, 0))
                item.setFont(font)
                item.setPos(float(loc), y_offset)
                self._plot.addItem(item)
                self._base_items.append(item)

        def _draw_trim_regions(self):
            """Shade trimmed (low-quality) regions."""
            if not self._chrom or len(self._chrom.peak_locations) == 0:
                return

            trim_color = "#ff444422" if self._dark_mode else "#ff000011"

            # Left trim region
            if self._chrom.trim_start > 0:
                left_end = int(self._chrom.peak_locations[self._chrom.trim_start])
                self._trim_left = LinearRegionItem(
                    values=[0, left_end],
                    brush=pg.mkBrush(trim_color),
                    pen=pg.mkPen(None),
                    movable=False,
                )
                self._plot.addItem(self._trim_left)

            # Right trim region
            if self._chrom.trim_end < self._chrom.num_bases:
                right_start = int(self._chrom.peak_locations[self._chrom.trim_end - 1])
                right_end = self._chrom.trace_length
                self._trim_right = LinearRegionItem(
                    values=[right_start, right_end],
                    brush=pg.mkBrush(trim_color),
                    pen=pg.mkPen(None),
                    movable=False,
                )
                self._plot.addItem(self._trim_right)

        def _auto_range(self):
            """Fit the view to show all data."""
            self._plot.autoRange()

        def scroll_to_base(self, base_index: int):
            """Center the view on a specific base position."""
            if not self._chrom or base_index >= len(self._chrom.peak_locations):
                return
            loc = int(self._chrom.peak_locations[base_index])
            vr = self._plot.getViewBox().viewRange()
            x_span = vr[0][1] - vr[0][0]
            self._plot.setXRange(loc - x_span / 2, loc + x_span / 2, padding=0)

        def highlight_bases(self, positions: list[int], color: str = "#ff0000"):
            """Highlight specific base positions (e.g., mismatches)."""
            for item in self._highlight_items:
                self._plot.removeItem(item)
            self._highlight_items.clear()

            if not self._chrom:
                return

            for pos in positions:
                if pos >= len(self._chrom.peak_locations):
                    continue
                loc = int(self._chrom.peak_locations[pos])
                line = InfiniteLine(
                    pos=loc, angle=90,
                    pen=pg.mkPen(color, width=2, style=Qt.PenStyle.DashLine),
                )
                self._plot.addItem(line)
                self._highlight_items.append(line)

        def set_dark_mode(self, dark: bool):
            """Switch between dark and light themes."""
            self._dark_mode = dark
            bg = "#1e1e2e" if dark else "#fafafa"
            self._plot_widget.setBackground(bg)
            if self._chrom:
                self.set_chromatogram(self._chrom)

        def _on_click(self, event):
            """Handle mouse click on the plot."""
            if not self._chrom or len(self._chrom.peak_locations) == 0:
                return

            pos = event.scenePos()
            mouse_point = self._plot.getViewBox().mapSceneToView(pos)
            x = mouse_point.x()

            # Find nearest base
            distances = np.abs(self._chrom.peak_locations.astype(float) - x)
            nearest = int(np.argmin(distances))

            if event.double():
                self._edit_base(nearest)
            else:
                self.base_selected.emit(nearest)
                # Move cursor
                loc = float(self._chrom.peak_locations[nearest])
                self._cursor_line.setValue(loc)

        def _edit_base(self, position: int):
            """Open inline editor for a base call."""
            if not self._chrom:
                return
            current = self._chrom.basecalls[position]
            new_base, ok = QInputDialog.getText(
                self, "Edit Base",
                f"Position {position + 1} (current: {current}):",
                text=current,
            )
            if ok and new_base:
                new_base = new_base.upper().strip()
                if len(new_base) == 1 and new_base in "ACGTNRYSWKMBDHV-":
                    self._chrom.edit_base(position, new_base)
                    self.base_edited.emit(position, new_base)
                    self.set_chromatogram(self._chrom)  # Redraw

        def _on_range_changed(self, viewBox, ranges):
            """Emit signal when view range changes (for syncing)."""
            x_range = ranges[0]
            self.view_range_changed.emit(x_range[0], x_range[1])

        def export_to_image(self, filepath: str):
            """Export the current view as an image."""
            from pyqtgraph.exporters import ImageExporter, SVGExporter
            if filepath.lower().endswith(".svg"):
                exporter = SVGExporter(self._plot)
            else:
                exporter = ImageExporter(self._plot)
                exporter.parameters()["width"] = 3200
            exporter.export(filepath)

else:
    # Stub for environments without PyQt6
    class TraceWidget:
        def __init__(self, *args, **kwargs):
            raise ImportError("PyQt6 and pyqtgraph required for TraceWidget")
