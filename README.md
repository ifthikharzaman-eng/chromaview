# ChromaView — Sanger Sequencing Chromatogram Analyzer

**ChromaView** is a modern, open-source desktop application for visualizing, editing, and analyzing Sanger DNA sequencing chromatogram files. It provides functionality comparable to Chromas, FinchTV, and SnapGene's trace viewer, built with a clean Python stack.

---

## Features

- **File Support**: Import `.ab1` and `.scf` chromatogram files
- **Chromatogram Visualization**: Interactive fluorescence traces (A/T/G/C) with zoom, pan, and scroll
- **Base Calling Display**: Called bases aligned to peaks with quality-score coloring
- **Manual Editing**: Click-to-edit base calls with undo/redo
- **Quality Analysis**: Phred score bar chart, signal-to-noise ratio, ambiguous base detection
- **Trimming**: Automatic and manual trimming of low-quality ends (modified Mott algorithm)
- **Sequence Operations**: Reverse complement, search by position or motif
- **Reference Comparison**: Align chromatogram against a reference sequence, highlight mismatches/SNPs
- **Batch Processing**: Process multiple files for quality reports
- **Export**: FASTA, CSV, publication-quality PNG/SVG chromatogram images
- **Modern UI**: Dark/light themes, dockable panels, synchronized sequence viewer

## Architecture

```
chromaview/
├── core/           # Data models and file parsers
│   ├── ab1_parser.py       # AB1/ABIF file parser
│   ├── scf_parser.py       # SCF file parser
│   ├── models.py           # Chromatogram data model
│   └── sequence_ops.py     # Reverse complement, alignment, search
├── gui/            # PyQt6 user interface
│   ├── main_window.py      # Main application window
│   ├── trace_widget.py     # Chromatogram trace canvas
│   ├── sequence_bar.py     # Synchronized sequence viewer
│   ├── quality_widget.py   # Phred score bar chart
│   ├── file_browser.py     # Sidebar file explorer
│   ├── dialogs.py          # Trim, compare, batch dialogs
│   └── theme.py            # Dark/light theme manager
├── analysis/       # Analysis algorithms
│   ├── quality.py          # Quality filtering, trimming
│   ├── peaks.py            # Peak detection and metrics
│   └── mutation.py         # SNP detection, mismatch calling
├── export/         # Output generators
│   ├── fasta.py            # FASTA export
│   ├── csv_export.py       # CSV export
│   └── image.py            # Publication-quality image export
├── utils/          # Shared utilities
│   └── helpers.py
├── tests/          # Unit tests
│   ├── test_ab1_parser.py
│   ├── test_models.py
│   ├── test_quality.py
│   └── test_sequence_ops.py
└── app.py          # Application entry point
```

## Installation

### Requirements

- Python 3.10+
- System: Windows, macOS, or Linux

### Install from source

```bash
git clone https://github.com/yourname/chromaview.git
cd chromaview
pip install -e .
```

### Run

```bash
chromaview
# or
python -m chromaview
```

### Install dependencies only

```bash
pip install -r requirements.txt
```

## Usage

1. **Open a file**: File → Open or drag-and-drop an `.ab1` file
2. **Navigate**: Scroll/zoom the chromatogram; click a base to select it
3. **Edit**: Double-click a base to change its call
4. **Trim**: Analysis → Trim Low-Quality Ends
5. **Compare**: Analysis → Compare to Reference
6. **Export**: File → Export as FASTA / CSV / Image

## Development

```bash
# Run tests
pytest chromaview/tests/ -v

# Run with debug logging
chromaview --debug
```

## Tech Stack

| Component       | Library        |
|----------------|----------------|
| GUI Framework  | PyQt6          |
| Chromatogram   | PyQtGraph      |
| AB1 Parsing    | Biopython      |
| Data handling  | NumPy          |
| Testing        | pytest         |

## License

MIT License. See [LICENSE](LICENSE).

## Changelog

See [CHANGELOG.md](CHANGELOG.md).
