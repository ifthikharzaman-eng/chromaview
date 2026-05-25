"""AB1 (Applied Biosystems) chromatogram file parser.

AB1 files use the ABIF (Applied Biosystems Inc. Format) binary format.
Structure:
    - 128-byte header containing magic bytes "ABIF" and a pointer to the
      tag directory.
    - Tag directory: array of 28-byte entries, each describing one data
      element by a 4-char tag name + tag number.

Key tags used:
    DATA 9-12  : Analyzed fluorescence traces (baseline-corrected)
    DATA 1-4   : Raw fluorescence traces
    FWO_ 1     : Filter Wheel Order — 4-char string mapping channel order
                 to nucleotides (e.g. "GATC" means channel 1=G, 2=A, 3=T, 4=C)
    PBAS 1     : Primary base calls (raw)
    PBAS 2     : Edited base calls
    PLOC 1     : Peak locations (sample indices for each called base)
    PLOC 2     : Edited peak locations
    PCON 1/2   : Quality confidence values (Phred scores)
    SMPL 1     : Sample name
    MODL 1     : Instrument model
    RUND 1/2   : Run date/time
    SPAC 1     : Mean peak spacing

We use Biopython's AbiIterator for robust ABIF parsing, falling back
to our own struct-based reader for tags Biopython doesn't expose.
"""
from __future__ import annotations

import struct
from pathlib import Path
from typing import Union

import numpy as np

from .models import Chromatogram, ChromatogramMetadata


def parse_ab1(filepath: Union[str, Path]) -> Chromatogram:
    """Parse an AB1 file and return a Chromatogram instance.

    Uses direct binary parsing of the ABIF directory for maximum control
    over which tags we extract. Biopython's SeqIO.read(handle, 'abi')
    is an alternative but doesn't expose raw traces or all metadata.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"AB1 file not found: {filepath}")

    with open(filepath, "rb") as fh:
        data = fh.read()

    # ── Validate header ──────────────────────────────────────────
    magic = data[:4]
    if magic != b"ABIF":
        raise ValueError(f"Not a valid AB1/ABIF file (magic={magic!r})")

    version = struct.unpack(">H", data[4:6])[0]

    # Directory pointer is at offset 26 in the header
    dir_tag = data[6:10]  # should be 'tdir'
    dir_entry_count = struct.unpack(">i", data[18:22])[0]
    dir_offset = struct.unpack(">i", data[26:30])[0]

    # ── Read directory entries ───────────────────────────────────
    tags = {}
    for i in range(dir_entry_count):
        entry_offset = dir_offset + i * 28
        entry = data[entry_offset : entry_offset + 28]
        tag_name = entry[0:4].decode("ascii", errors="replace")
        tag_number = struct.unpack(">i", entry[4:8])[0]
        element_type = struct.unpack(">H", entry[8:10])[0]
        element_size = struct.unpack(">H", entry[10:12])[0]
        num_elements = struct.unpack(">i", entry[12:16])[0]
        data_size = struct.unpack(">i", entry[16:20])[0]

        # Data location: if total size ≤ 4 bytes, data is stored inline
        # at offset 20; otherwise offset 20 holds a pointer
        if data_size <= 4:
            tag_data = entry[20 : 20 + data_size]
        else:
            data_offset = struct.unpack(">i", entry[20:24])[0]
            tag_data = data[data_offset : data_offset + data_size]

        key = (tag_name, tag_number)
        tags[key] = {
            "type": element_type,
            "size": element_size,
            "count": num_elements,
            "data": tag_data,
        }

    # ── Helper: extract array of shorts ──────────────────────────
    def get_short_array(tag_name: str, tag_number: int) -> np.ndarray:
        key = (tag_name, tag_number)
        if key not in tags:
            return np.array([], dtype=np.int16)
        raw = tags[key]["data"]
        count = tags[key]["count"]
        arr = np.array(struct.unpack(f">{count}h", raw[:count * 2]), dtype=np.int16)
        return arr

    def get_string(tag_name: str, tag_number: int) -> str:
        key = (tag_name, tag_number)
        if key not in tags:
            return ""
        raw = tags[key]["data"]
        # Strings are typically Pascal-style (length-prefixed) or null-terminated
        try:
            text = raw.decode("ascii", errors="replace").strip("\x00").strip()
        except Exception:
            text = ""
        return text

    def get_float(tag_name: str, tag_number: int) -> float:
        key = (tag_name, tag_number)
        if key not in tags:
            return 0.0
        raw = tags[key]["data"]
        if len(raw) >= 4:
            return struct.unpack(">f", raw[:4])[0]
        return 0.0

    # ── Extract channel order ────────────────────────────────────
    fwo = get_string("FWO_", 1)
    if len(fwo) < 4:
        fwo = "GATC"  # common default
    channel_order = list(fwo[:4])  # e.g. ['G','A','T','C']

    # ── Extract analyzed traces (DATA 9-12) ──────────────────────
    traces = {}
    for i, base in enumerate(channel_order):
        arr = get_short_array("DATA", 9 + i)
        traces[base.upper()] = arr.astype(np.float64)

    # ── Extract raw traces (DATA 1-4) ────────────────────────────
    raw_traces = {}
    for i, base in enumerate(channel_order):
        arr = get_short_array("DATA", 1 + i)
        raw_traces[base.upper()] = arr.astype(np.float64)

    # ── Base calls ───────────────────────────────────────────────
    # Prefer edited calls (PBAS 2), fall back to raw (PBAS 1)
    basecall_str = get_string("PBAS", 2) or get_string("PBAS", 1)
    basecalls = list(basecall_str.upper())

    # ── Peak locations ───────────────────────────────────────────
    # PLOC 2 = edited, PLOC 1 = raw
    peak_locs = get_short_array("PLOC", 2)
    if len(peak_locs) == 0:
        peak_locs = get_short_array("PLOC", 1)
    peak_locations = peak_locs.astype(int)

    # ── Quality scores (Phred) ───────────────────────────────────
    # PCON stores 1-byte unsigned ints, not shorts
    quality = np.array([], dtype=int)
    for pcon_num in (2, 1):
        key = ("PCON", pcon_num)
        if key in tags:
            raw = tags[key]["data"]
            count = tags[key]["count"]
            quality = np.array(
                struct.unpack(f">{count}B", raw[:count]), dtype=int
            )
            break

    # Ensure arrays match in length
    n_bases = len(basecalls)
    if len(peak_locations) > n_bases:
        peak_locations = peak_locations[:n_bases]
    if len(quality) > n_bases:
        quality = quality[:n_bases]
    elif len(quality) < n_bases:
        quality = np.pad(quality, (0, n_bases - len(quality)), constant_values=0)

    # ── Metadata ─────────────────────────────────────────────────
    meta = ChromatogramMetadata(
        sample_name=get_string("SMPL", 1),
        file_path=str(filepath),
        file_format="ab1",
        instrument_model=get_string("MODL", 1),
        spacing=get_float("SPAC", 1),
        comment=get_string("CMNT", 1),
    )

    # ── Assemble Chromatogram ────────────────────────────────────
    chrom = Chromatogram(
        traces=traces,
        raw_traces=raw_traces,
        basecalls=basecalls,
        peak_locations=peak_locations,
        quality_scores=quality,
        metadata=meta,
    )
    return chrom
