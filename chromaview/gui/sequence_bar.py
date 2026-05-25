"""Sequence viewer bar — synchronized with the chromatogram trace.

Displays the called sequence as a horizontal strip of colored base
letters. Stays in sync with the trace widget's x-range. Edited bases
are underlined, trimmed regions are dimmed.
"""
from __future__ import annotations

from typing import Optional

try:
    from PyQt6.QtCore import Qt, pyqtSignal
    from PyQt6.QtGui import QFont, QFontMetrics, QPainter, QColor, QPen
    from PyQt6.QtWidgets import QWidget, QHBoxLayout, QScrollBar
    HAS_GUI = True
except ImportError:
    HAS_GUI = False

from ..core.models import Chromatogram, BASE_COLORS, BASE_COLORS_DARK


if HAS_GUI:

    class SequenceBar(QWidget):
        """Horizontal sequence display that syncs with TraceWidget.

        Signals:
            base_clicked(int): Emitted when a base letter is clicked.
        """

        base_clicked = pyqtSignal(int)

        def __init__(self, parent=None, dark_mode: bool = True):
            super().__init__(parent)
            self._chrom: Optional[Chromatogram] = None
            self._dark_mode = dark_mode
            self._view_start = 0  # first visible base index
            self._view_end = 50   # last visible base index
            self._selected_base = -1
            self._highlight_positions: set[int] = set()

            self._font = QFont("Consolas", 11, QFont.Weight.Bold)
            self._metrics = QFontMetrics(self._font)
            self._char_width = self._metrics.horizontalAdvance("A") + 4

            self.setMinimumHeight(36)
            self.setMaximumHeight(36)
            self.setCursor(Qt.CursorShape.PointingHandCursor)

        def set_chromatogram(self, chrom: Chromatogram):
            self._chrom = chrom
            self._view_start = 0
            self._view_end = min(chrom.num_bases, 100)
            self.update()

        def set_view_range(self, sample_start: float, sample_end: float):
            """Update visible range based on trace widget's x-range (sample coords)."""
            if not self._chrom or len(self._chrom.peak_locations) == 0:
                return

            locs = self._chrom.peak_locations.astype(float)
            # Find base indices within the sample range
            mask = (locs >= sample_start) & (locs <= sample_end)
            visible = mask.nonzero()[0]
            if len(visible) > 0:
                self._view_start = int(visible[0])
                self._view_end = int(visible[-1]) + 1
            self.update()

        def set_selected_base(self, index: int):
            self._selected_base = index
            self.update()

        def set_highlights(self, positions: list[int]):
            """Highlight specific positions (mismatches, search hits)."""
            self._highlight_positions = set(positions)
            self.update()

        def set_dark_mode(self, dark: bool):
            self._dark_mode = dark
            self.update()

        def paintEvent(self, event):
            if not self._chrom:
                return

            painter = QPainter(self)
            painter.setRenderHint(QPainter.RenderHint.Antialiasing)
            painter.setFont(self._font)

            colors = BASE_COLORS_DARK if self._dark_mode else BASE_COLORS
            bg = QColor("#2a2a3d") if self._dark_mode else QColor("#ffffff")
            text_dim = QColor("#555555") if self._dark_mode else QColor("#bbbbbb")

            painter.fillRect(self.rect(), bg)

            # Calculate layout
            available_width = self.width()
            n_visible = self._view_end - self._view_start
            if n_visible <= 0:
                painter.end()
                return

            char_w = max(available_width / n_visible, 8)
            y_center = self.height() // 2 + 5

            for i in range(self._view_start, min(self._view_end, self._chrom.num_bases)):
                x = (i - self._view_start) * char_w
                base = self._chrom.basecalls[i]
                is_trimmed = i < self._chrom.trim_start or i >= self._chrom.trim_end

                # Background for selected base
                if i == self._selected_base:
                    painter.fillRect(
                        int(x), 0, int(char_w), self.height(),
                        QColor("#7c9ff544" if self._dark_mode else "#4a6cf733"),
                    )

                # Background for highlighted bases
                if i in self._highlight_positions:
                    painter.fillRect(
                        int(x), 0, int(char_w), self.height(),
                        QColor("#ff444444"),
                    )

                # Base letter color
                if is_trimmed:
                    color = text_dim
                else:
                    color = QColor(colors.get(base.upper(), "#888888"))

                painter.setPen(color)
                painter.drawText(
                    int(x + char_w / 2 - self._metrics.horizontalAdvance(base) / 2),
                    y_center, base,
                )

                # Position number every 10 bases
                if (i + 1) % 10 == 0:
                    painter.setPen(QColor("#666666"))
                    num_str = str(i + 1)
                    painter.setFont(QFont("Consolas", 7))
                    painter.drawText(
                        int(x + char_w / 2 - self._metrics.horizontalAdvance(num_str) / 4),
                        10, num_str,
                    )
                    painter.setFont(self._font)

                # Edited base marker (dot below)
                if (self._chrom._original_basecalls and
                        i < len(self._chrom._original_basecalls) and
                        self._chrom._original_basecalls[i] != self._chrom.basecalls[i]):
                    painter.setPen(QPen(QColor("#ff8800"), 2))
                    painter.drawEllipse(int(x + char_w / 2 - 2), self.height() - 6, 4, 4)

            painter.end()

        def mousePressEvent(self, event):
            if not self._chrom:
                return
            n_visible = self._view_end - self._view_start
            if n_visible <= 0:
                return
            char_w = self.width() / n_visible
            index = self._view_start + int(event.position().x() / char_w)
            if 0 <= index < self._chrom.num_bases:
                self._selected_base = index
                self.base_clicked.emit(index)
                self.update()

else:
    class SequenceBar:
        def __init__(self, *args, **kwargs):
            raise ImportError("PyQt6 required for SequenceBar")
