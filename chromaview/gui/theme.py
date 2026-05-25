"""Theme manager for dark/light mode switching.

Provides a QSS stylesheet and color constants that all widgets
reference for consistent appearance. Scientific software convention:
dark background for trace viewing (reduces eye strain during analysis),
light mode for printing/export.
"""
from __future__ import annotations

from dataclasses import dataclass

from ..core.models import BASE_COLORS, BASE_COLORS_DARK


@dataclass
class ThemeColors:
    """Color palette for the current theme."""

    background: str
    surface: str
    text: str
    text_secondary: str
    border: str
    accent: str
    quality_good: str
    quality_ok: str
    quality_bad: str
    selection: str
    trim_region: str
    base_colors: dict[str, str]


DARK_THEME = ThemeColors(
    background="#1e1e2e",
    surface="#2a2a3d",
    text="#e0e0e0",
    text_secondary="#8888aa",
    border="#3a3a50",
    accent="#7c9ff5",
    quality_good="#2ecc71",
    quality_ok="#f39c12",
    quality_bad="#e74c3c",
    selection="#44447788",
    trim_region="#ff444422",
    base_colors=BASE_COLORS_DARK,
)

LIGHT_THEME = ThemeColors(
    background="#fafafa",
    surface="#ffffff",
    text="#1a1a2e",
    text_secondary="#666688",
    border="#d0d0dd",
    accent="#4a6cf7",
    quality_good="#27ae60",
    quality_ok="#e67e22",
    quality_bad="#c0392b",
    selection="#4a6cf733",
    trim_region="#ff000011",
    base_colors=BASE_COLORS,
)


def get_stylesheet(theme: ThemeColors) -> str:
    """Generate QSS stylesheet for the given theme."""
    return f"""
    QMainWindow {{
        background-color: {theme.background};
    }}
    QWidget {{
        background-color: {theme.background};
        color: {theme.text};
        font-family: "Segoe UI", "SF Pro Text", "Noto Sans", sans-serif;
        font-size: 13px;
    }}
    QDockWidget {{
        titlebar-close-icon: none;
        color: {theme.text};
        font-weight: bold;
        border: 1px solid {theme.border};
    }}
    QDockWidget::title {{
        background-color: {theme.surface};
        padding: 6px;
        border-bottom: 1px solid {theme.border};
    }}
    QMenuBar {{
        background-color: {theme.surface};
        color: {theme.text};
        border-bottom: 1px solid {theme.border};
        padding: 2px;
    }}
    QMenuBar::item:selected {{
        background-color: {theme.accent};
        color: white;
        border-radius: 4px;
    }}
    QMenu {{
        background-color: {theme.surface};
        color: {theme.text};
        border: 1px solid {theme.border};
        padding: 4px;
    }}
    QMenu::item:selected {{
        background-color: {theme.accent};
        color: white;
    }}
    QToolBar {{
        background-color: {theme.surface};
        border-bottom: 1px solid {theme.border};
        padding: 3px;
        spacing: 4px;
    }}
    QPushButton {{
        background-color: {theme.surface};
        color: {theme.text};
        border: 1px solid {theme.border};
        border-radius: 4px;
        padding: 5px 14px;
        font-weight: 500;
    }}
    QPushButton:hover {{
        background-color: {theme.accent};
        color: white;
        border-color: {theme.accent};
    }}
    QPushButton:pressed {{
        background-color: {theme.accent}cc;
    }}
    QLineEdit {{
        background-color: {theme.surface};
        color: {theme.text};
        border: 1px solid {theme.border};
        border-radius: 4px;
        padding: 5px 8px;
    }}
    QLineEdit:focus {{
        border-color: {theme.accent};
    }}
    QTreeView, QListView {{
        background-color: {theme.surface};
        color: {theme.text};
        border: 1px solid {theme.border};
        alternate-background-color: {theme.background};
        outline: none;
    }}
    QTreeView::item:selected, QListView::item:selected {{
        background-color: {theme.accent};
        color: white;
    }}
    QHeaderView::section {{
        background-color: {theme.surface};
        color: {theme.text};
        border: none;
        border-right: 1px solid {theme.border};
        border-bottom: 1px solid {theme.border};
        padding: 4px 8px;
        font-weight: bold;
    }}
    QStatusBar {{
        background-color: {theme.surface};
        color: {theme.text_secondary};
        border-top: 1px solid {theme.border};
    }}
    QScrollBar:horizontal {{
        height: 10px;
        background: {theme.background};
    }}
    QScrollBar::handle:horizontal {{
        background: {theme.border};
        border-radius: 5px;
        min-width: 30px;
    }}
    QScrollBar:vertical {{
        width: 10px;
        background: {theme.background};
    }}
    QScrollBar::handle:vertical {{
        background: {theme.border};
        border-radius: 5px;
        min-height: 30px;
    }}
    QTabWidget::pane {{
        border: 1px solid {theme.border};
    }}
    QTabBar::tab {{
        background-color: {theme.surface};
        color: {theme.text_secondary};
        border: 1px solid {theme.border};
        border-bottom: none;
        padding: 6px 16px;
        margin-right: 2px;
    }}
    QTabBar::tab:selected {{
        background-color: {theme.background};
        color: {theme.accent};
        font-weight: bold;
    }}
    QSplitter::handle {{
        background-color: {theme.border};
    }}
    QLabel#status_label {{
        color: {theme.text_secondary};
        font-size: 12px;
    }}
    QDialog {{
        background-color: {theme.background};
    }}
    QGroupBox {{
        border: 1px solid {theme.border};
        border-radius: 4px;
        margin-top: 8px;
        padding-top: 16px;
        font-weight: bold;
    }}
    QGroupBox::title {{
        subcontrol-origin: margin;
        left: 10px;
        padding: 0 4px;
    }}
    QSpinBox, QDoubleSpinBox {{
        background-color: {theme.surface};
        color: {theme.text};
        border: 1px solid {theme.border};
        border-radius: 4px;
        padding: 3px;
    }}
    QTextEdit {{
        background-color: {theme.surface};
        color: {theme.text};
        border: 1px solid {theme.border};
        border-radius: 4px;
        font-family: "JetBrains Mono", "Fira Code", "Consolas", monospace;
    }}
    """
