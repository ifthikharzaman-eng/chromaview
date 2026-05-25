"""SCF (Standard Chromatogram Format) file parser.

SCF is a simpler, standardized format for chromatogram data.
Structure (v3):
    - 128-byte header with counts, offsets, and version
    - Sample data section: interleaved or sequential trace points
    - Base data section: peak indices, accuracies, and called bases
    - Comment section: key=value pairs

SCF v3 stores samples as delta-delta encoded uint8 or uint16 values
for compression. We decode these back to absolute intensities.
"""
from __future__ import annotations

import struct
from pathlib import Path
from typing import Union

import numpy as np

from .models import Chromatogram, ChromatogramMetadata


def _decode_delta_delta(data: np.ndarray) -> np.ndarray:
    """Reverse SCF v3 delta-delta encoding.

    SCF v3 stores trace values as second-order differences.
    To reconstruct: first cumulative sum reverses the outer delta,
    second cumulative sum reverses the inner delta.
    """
    result = np.cumsum(data).astype(np.float64)
    result = np.cumsum(result).astype(np.float64)
    return result


def parse_scf(filepath: Union[str, Path]) -> Chromatogram:
    """Parse an SCF file and return a Chromatogram instance."""
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"SCF file not found: {filepath}")

    with open(filepath, "rb") as fh:
        data = fh.read()

    # ── Header (128 bytes) ───────────────────────────────────────
    magic = data[0:4]
    if magic not in (b".scf", b"SCF\x03"):
        # Some SCF files have different magic; try to proceed
        pass

    (
        num_samples,
        samples_offset,
        num_bases,
        _bases_left_clip,
        _bases_right_clip,
        bases_offset,
        comments_size,
        comments_offset,
        version_str,
        sample_size,
    ) = struct.unpack(">IIIIIIIi4sI", data[4:44])

    # Decode version
    try:
        version = float(version_str.decode("ascii", errors="replace").strip("\x00"))
    except (ValueError, UnicodeDecodeError):
        version = 3.0

    # ── Read trace samples ───────────────────────────────────────
    # Channel order in SCF is always A, C, G, T
    channel_map = ["A", "C", "G", "T"]
    traces = {}

    if version >= 3.0:
        # v3: samples stored sequentially by channel, then delta-delta encoded
        if sample_size == 1:
            fmt_char = "B"
            fmt_size = 1
        else:
            fmt_char = ">H"
            fmt_size = 2

        for ch_idx, base in enumerate(channel_map):
            offset = samples_offset + ch_idx * num_samples * fmt_size
            if fmt_size == 1:
                raw = np.array(
                    struct.unpack(f">{num_samples}B",
                                 data[offset:offset + num_samples]),
                    dtype=np.int16,
                )
            else:
                raw = np.array(
                    struct.unpack(f">{num_samples}H",
                                 data[offset:offset + num_samples * 2]),
                    dtype=np.int16,
                )
            traces[base] = _decode_delta_delta(raw)
    else:
        # v2: interleaved samples (A1,C1,G1,T1, A2,C2,G2,T2, ...)
        all_samples = {b: [] for b in channel_map}
        for i in range(num_samples):
            offset = samples_offset + i * 8  # 4 channels × 2 bytes
            for ch_idx, base in enumerate(channel_map):
                val = struct.unpack(">H", data[offset + ch_idx * 2:offset + ch_idx * 2 + 2])[0]
                all_samples[base].append(val)
        for base in channel_map:
            traces[base] = np.array(all_samples[base], dtype=np.float64)

    # ── Read base calls ──────────────────────────────────────────
    basecalls = []
    peak_locations = []
    quality_a = []
    quality_c = []
    quality_g = []
    quality_t = []

    if version >= 3.0:
        # v3 base records: 12 bytes each
        # peak_index(4) + prob_A(1) + prob_C(1) + prob_G(1) + prob_T(1)
        # + base(1) + spare(3)
        # But typically laid out as sequential arrays
        peaks_offset = bases_offset
        peaks_data = data[peaks_offset:peaks_offset + num_bases * 4]
        peak_locations = np.array(
            struct.unpack(f">{num_bases}I", peaks_data), dtype=int
        )

        # Accuracy arrays (1 byte each, 4 channels)
        acc_offset = peaks_offset + num_bases * 4
        prob_a = struct.unpack(f">{num_bases}B", data[acc_offset:acc_offset + num_bases])
        acc_offset += num_bases
        prob_c = struct.unpack(f">{num_bases}B", data[acc_offset:acc_offset + num_bases])
        acc_offset += num_bases
        prob_g = struct.unpack(f">{num_bases}B", data[acc_offset:acc_offset + num_bases])
        acc_offset += num_bases
        prob_t = struct.unpack(f">{num_bases}B", data[acc_offset:acc_offset + num_bases])
        acc_offset += num_bases

        # Called bases
        bases_raw = data[acc_offset:acc_offset + num_bases]
        basecalls = list(bases_raw.decode("ascii", errors="replace").upper())

        # Quality = max of the 4 probabilities per position
        quality_scores = np.maximum(
            np.maximum(prob_a, prob_c),
            np.maximum(prob_g, prob_t),
        )
        quality_scores = np.array(quality_scores, dtype=int)
    else:
        # v2: 12-byte base records, interleaved
        quality_scores_list = []
        for i in range(num_bases):
            offset = bases_offset + i * 12
            peak_idx = struct.unpack(">I", data[offset:offset + 4])[0]
            peak_locations.append(peak_idx)
            pa, pc, pg, pt = struct.unpack("4B", data[offset + 4:offset + 8])
            base_char = chr(data[offset + 8])
            basecalls.append(base_char.upper())
            quality_scores_list.append(max(pa, pc, pg, pt))
        peak_locations = np.array(peak_locations, dtype=int)
        quality_scores = np.array(quality_scores_list, dtype=int)

    # ── Comments ─────────────────────────────────────────────────
    comment_text = ""
    if comments_size > 0:
        raw_comments = data[comments_offset:comments_offset + comments_size]
        try:
            comment_text = raw_comments.decode("ascii", errors="replace").strip("\x00")
        except Exception:
            pass

    # ── Metadata ─────────────────────────────────────────────────
    meta = ChromatogramMetadata(
        sample_name=filepath.stem,
        file_path=str(filepath),
        file_format="scf",
        comment=comment_text,
    )

    return Chromatogram(
        traces=traces,
        raw_traces={},
        basecalls=basecalls,
        peak_locations=peak_locations,
        quality_scores=quality_scores,
        metadata=meta,
    )
