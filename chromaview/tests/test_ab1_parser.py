"""Tests for the AB1 parser.

Since we can't ship a real AB1 file in the test suite, we test:
1. Header validation (reject non-ABIF files)
2. Synthetic AB1 construction and parsing
"""
import struct
import tempfile
from pathlib import Path

import numpy as np
import pytest

from chromaview.core.ab1_parser import parse_ab1


class TestAB1Parser:
    def test_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            parse_ab1("/nonexistent/file.ab1")

    def test_invalid_magic(self):
        with tempfile.NamedTemporaryFile(suffix=".ab1", delete=False) as f:
            f.write(b"NOT_ABIF" + b"\x00" * 200)
            f.flush()
            with pytest.raises(ValueError, match="Not a valid AB1"):
                parse_ab1(f.name)

    def test_minimal_synthetic_ab1(self):
        """Build a minimal valid ABIF file and parse it."""
        # This creates a bare-minimum ABIF with:
        # - Valid header
        # - FWO_ tag (channel order)
        # - DATA 9-12 (analyzed traces, 100 points each)
        # - PBAS 1 (base calls)
        # - PLOC 1 (peak locations)
        # - PCON 1 (quality scores)
        # - SMPL 1 (sample name)

        n_samples = 100
        n_bases = 10
        bases = b"ATGCATGCAT"
        peak_locs = np.linspace(10, 90, n_bases).astype(np.int16)
        quality = np.full(n_bases, 30, dtype=np.uint8)
        trace = np.random.default_rng(42).integers(0, 500, n_samples).astype(np.int16)

        # Build tag directory entries
        tags = []
        data_blobs = []
        data_offset_base = 256 + 28 * 9  # header(128) + padding(128) + 9 entries × 28

        def add_tag(name, number, element_type, element_size, count, raw_data):
            tag_data = raw_data
            data_size = len(tag_data)
            entry = struct.pack(">4siHHiI",
                                name.encode("ascii"), number,
                                element_type, element_size,
                                count, data_size)
            # Data offset or inline
            if data_size <= 4:
                entry += tag_data.ljust(4, b"\x00")
            else:
                offset = data_offset_base + sum(len(d) for d in data_blobs)
                entry += struct.pack(">I", offset)
                data_blobs.append(tag_data)
            entry += b"\x00" * (28 - len(entry))  # pad to 28 bytes
            tags.append(entry)

        # FWO_ 1: channel order "GATC" (4 bytes, stored inline)
        add_tag("FWO_", 1, 2, 1, 4, b"GATC")

        # DATA 9-12: traces (4 channels)
        for ch in range(4):
            raw = struct.pack(f">{n_samples}h", *trace)
            add_tag("DATA", 9 + ch, 4, 2, n_samples, raw)

        # PBAS 1: base calls
        add_tag("PBAS", 1, 2, 1, n_bases, bases)

        # PLOC 1: peak locations
        ploc_raw = struct.pack(f">{n_bases}h", *peak_locs)
        add_tag("PLOC", 1, 4, 2, n_bases, ploc_raw)

        # PCON 1: quality scores
        pcon_raw = struct.pack(f">{n_bases}B", *quality)
        add_tag("PCON", 1, 1, 1, n_bases, pcon_raw)

        # SMPL 1: sample name
        sample_name = b"TestSample"
        add_tag("SMPL", 1, 18, 1, len(sample_name), sample_name)

        # Build header (128 bytes)
        n_entries = len(tags)
        dir_offset = 256  # after header + padding
        header = b"ABIF"
        header += struct.pack(">H", 101)  # version
        header += b"tdir"  # dir tag name
        header += struct.pack(">i", 1)  # dir tag number
        header += struct.pack(">HH", 1023, 28)  # element type, size
        header += struct.pack(">i", n_entries)
        header += struct.pack(">I", n_entries * 28)  # data size
        header += struct.pack(">i", dir_offset)
        header = header.ljust(128, b"\x00")

        # Assemble file
        padding = b"\x00" * 128
        dir_data = b"".join(tags)
        blob_data = b"".join(data_blobs)
        file_data = header + padding + dir_data + blob_data

        with tempfile.NamedTemporaryFile(suffix=".ab1", delete=False) as f:
            f.write(file_data)
            f.flush()
            chrom = parse_ab1(f.name)

        assert chrom.num_bases == n_bases
        assert chrom.metadata.sample_name == "TestSample"
        assert len(chrom.traces) == 4
        assert "G" in chrom.traces
        assert "A" in chrom.traces
        assert "T" in chrom.traces
        assert "C" in chrom.traces
        assert chrom.sequence == "ATGCATGCAT"
        assert len(chrom.quality_scores) == n_bases
        assert all(q == 30 for q in chrom.quality_scores)
