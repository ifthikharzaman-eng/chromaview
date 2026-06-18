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
- **NCBI BLAST**: In-application `blastn` search against selectable nucleotide databases (`core_nt`, `nt`, `refseq_rna`), running on a background thread so the UI stays responsive. Results are displayed in a sortable table (accession, description, query coverage, identity, E-value, bit score) with a pairwise alignment view and direct links to NCBI records. A contact email for NCBI can be saved in application settings.
- **Forward + Reverse Consensus**: Quality-weighted consensus from a forward/reverse read pair. The reverse read is automatically reverse-complemented; the overlap is found by local alignment. Each consensus position is called from the higher-quality base, or an IUPAC ambiguity code when quality is equal. Assembly is gated on a minimum overlap identity threshold, and a clear diagnostic message is shown when the reads do not overlap.

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
│   ├── blast_dialog.py     # NCBI BLAST search dialog
│   ├── consensus_dialog.py # Forward+reverse consensus dialog
│   └── theme.py            # Dark/light theme manager
├── analysis/       # Analysis algorithms
│   ├── quality.py          # Quality filtering, trimming
│   ├── peaks.py            # Peak detection and metrics
│   ├── mutation.py         # SNP detection, mismatch calling
│   ├── blast.py            # NCBI BLAST integration (background thread)
│   └── consensus.py        # Forward+reverse quality-weighted consensus
├── export/         # Output generators
│   ├── fasta.py            # FASTA export
│   ├── csv_export.py       # CSV export
│   └── image.py            # Publication-quality image export
├── utils/          # Shared utilities
├── tests/          # Unit tests
│   ├── test_ab1_parser.py
│   ├── test_models.py
│   ├── test_quality.py
│   ├── test_sequence_ops.py
│   ├── test_blast.py
│   └── test_consensus.py
└── app.py          # Application entry point
```

## Installation

### Requirements

- Python 3.10+
- System: Windows, macOS, or Linux

### Install from source

```bash
git clone https://github.com/ifthikharzaman-eng/chromaview.git
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

### Standalone Windows executable

A self-contained Windows `.exe` can be built with PyInstaller using the provided helper script:

```bat
build_windows.bat
```

The script regenerates the application icon, runs PyInstaller against `ChromaView.spec`, and writes the output to `dist\ChromaView\ChromaView.exe`.

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
