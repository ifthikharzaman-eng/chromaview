"""ChromaView application entry point.

Usage:
    chromaview              # Launch GUI
    chromaview --debug      # Launch with debug logging
    python -m chromaview    # Alternative launch
"""
from __future__ import annotations

import sys
import logging
from pathlib import Path


def _resource_path(relative: str) -> Path:
    """Resolve a resource path for both dev and PyInstaller frozen builds."""
    if getattr(sys, "frozen", False) and hasattr(sys, "_MEIPASS"):
        return Path(sys._MEIPASS) / relative
    # Dev: resources/ lives next to the chromaview package directory
    return Path(__file__).parent.parent / relative


def main():
    """Launch the ChromaView application."""
    # Parse args
    debug = "--debug" in sys.argv

    logging.basicConfig(
        level=logging.DEBUG if debug else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )
    logger = logging.getLogger("chromaview")

    try:
        from PyQt6.QtWidgets import QApplication
        from PyQt6.QtCore import Qt
    except ImportError:
        print(
            "ERROR: PyQt6 is required but not installed.\n"
            "Install it with: pip install PyQt6",
            file=sys.stderr,
        )
        sys.exit(1)

    try:
        import pyqtgraph  # noqa: F401
    except ImportError:
        print(
            "ERROR: pyqtgraph is required but not installed.\n"
            "Install it with: pip install pyqtgraph",
            file=sys.stderr,
        )
        sys.exit(1)

    app = QApplication(sys.argv)
    app.setApplicationName("ChromaView")
    app.setApplicationVersion("0.1.0")
    app.setOrganizationName("ChromaView")

    # Window / taskbar icon
    from PyQt6.QtGui import QIcon
    icon_path = _resource_path("resources/icon.ico")
    if icon_path.exists():
        app.setWindowIcon(QIcon(str(icon_path)))

    # Enable high-DPI scaling
    app.setStyle("Fusion")

    from .gui.main_window import MainWindow

    window = MainWindow()
    window.show()

    # If a file was passed as argument, open it
    for arg in sys.argv[1:]:
        if not arg.startswith("--") and (arg.endswith(".ab1") or arg.endswith(".scf")):
            window._load_file(arg)
            break

    logger.info("ChromaView started")
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
