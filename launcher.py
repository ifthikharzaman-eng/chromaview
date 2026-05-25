"""PyInstaller entry point for ChromaView.

This top-level script avoids relative-import issues when frozen.
"""
from chromaview.app import main

if __name__ == "__main__":
    main()
