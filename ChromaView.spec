# -*- mode: python ; coding: utf-8 -*-
# ChromaView PyInstaller spec  —  onedir, no console window
# Rebuild: pyinstaller --noconfirm ChromaView.spec
#          (or double-click build_windows.bat)

from PyInstaller.utils.hooks import (
    collect_all,
    collect_submodules,
    collect_data_files,
)

# ── Collect everything from packages with dynamic imports ─────────────────

bio_datas,     bio_bins,     bio_hidden     = collect_all("Bio")
pyqtg_datas,   pyqtg_bins,   pyqtg_hidden   = collect_all("pyqtgraph")
chroma_hidden  = collect_submodules("chromaview")

# PyQt6 platform plugin is handled by PyInstaller's built-in hook;
# add a few extras that static analysis can miss.
extra_hidden = [
    "PyQt6.QtCore",
    "PyQt6.QtGui",
    "PyQt6.QtWidgets",
    "PyQt6.QtNetwork",
    "PyQt6.QtOpenGL",
    "PyQt6.QtOpenGLWidgets",
    "PyQt6.sip",
    # pyqtgraph needs these Qt bindings at runtime
    "PyQt6.QtSvg",
    "PyQt6.QtPrintSupport",
    # biopython XML parser (C extension)
    "Bio.Blast._parser",
    "Bio.Align._parser",
    "Bio.Blast.NCBIWWW",
    "Bio.SeqRecord",
    "Bio.SeqIO",
    "Bio.Seq",
    # numpy / scipy internals missed by analysis
    "numpy",
    "numpy.core._multiarray_umath",
]

a = Analysis(
    ["launcher.py"],
    pathex=[],
    binaries=bio_bins + pyqtg_bins,
    datas=[
        ("resources/icon.ico", "resources"),
    ] + bio_datas + pyqtg_datas,
    hiddenimports=bio_hidden + pyqtg_hidden + chroma_hidden + extra_hidden,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[
        "tkinter",
        "matplotlib",
        "IPython",
        "jupyter",
        "PIL",          # Pillow only needed for icon generation, not runtime
        "pytest",
        "setuptools",
    ],
    noarchive=False,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="ChromaView",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,          # UPX can break Qt DLLs; keep off
    console=False,      # no black console window
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon="resources/icon.ico",
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name="ChromaView",
)
