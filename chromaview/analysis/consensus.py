"""Consensus assembly from paired forward + reverse Sanger reads.

Algorithm:
1. Reverse-complement the reverse read (sequence + quality array reversed).
2. Smith-Waterman local alignment (Bio.Align.PairwiseAligner) finds the
   best-scoring overlap between trimmed forward and RC-reverse sequences.
3. Per-position consensus in the overlap:
   - Agree: use the base; quality = min(60, q_fwd + q_rc).
   - Disagree, quality difference >= threshold: take the higher-quality base.
   - Disagree, comparable/low quality: emit IUPAC ambiguity code.
4. Prefix (fwd-only) and suffix (rc-only) are taken verbatim.
5. Warn if overlap is absent or reads look mis-paired.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from ..core.sequence_ops import reverse_complement


# ── IUPAC ambiguity codes for two-base disagreements ─────────────────────────

_IUPAC: dict[frozenset, str] = {
    frozenset("AC"): "M",
    frozenset("AG"): "R",
    frozenset("AT"): "W",
    frozenset("CG"): "S",
    frozenset("CT"): "Y",
    frozenset("GT"): "K",
    frozenset("ACG"): "V",
    frozenset("ACT"): "H",
    frozenset("AGT"): "D",
    frozenset("CGT"): "B",
    frozenset("ACGT"): "N",
}


def _iupac(b1: str, b2: str) -> str:
    """IUPAC ambiguity code for two disagreeing unambiguous bases."""
    return _IUPAC.get(frozenset([b1.upper(), b2.upper()]), "N")


# ── Data classes ─────────────────────────────────────────────────────────────

@dataclass
class OverlapResult:
    """Smith-Waterman overlap between the forward and RC-reverse sequences."""
    fwd_start: int       # 0-based start of overlap in fwd (absolute index)
    fwd_end: int         # 0-based exclusive end of overlap in fwd
    rc_start: int        # 0-based start of overlap in RC sequence
    rc_end: int          # 0-based exclusive end of overlap in RC sequence
    fwd_aligned: str     # fwd substring with gap chars in the alignment
    rc_aligned: str      # rc substring with gap chars in the alignment
    overlap_length: int  # number of alignment columns (includes gap cols)
    identity: float      # fraction of identical non-gap columns
    score: float         # Smith-Waterman score
    has_overlap: bool    # True when overlap meets minimum thresholds


@dataclass
class ConsensusPosition:
    """One position in the assembled consensus."""
    base: str       # nucleotide letter (may be IUPAC ambiguity code)
    quality: int    # Phred score
    source: str     # "fwd" | "rev" | "agree" | "disagree" | "ambiguous"


@dataclass
class ConsensusResult:
    """Full result of assembling two Sanger reads."""
    consensus_seq: str
    consensus_quality: list[int]
    positions: list[ConsensusPosition]
    overlap: OverlapResult
    fwd_only_length: int   # bases contributed only by forward (before overlap)
    rev_only_length: int   # bases contributed only by RC reverse (after overlap)
    disagreements: int     # overlap positions where one base won on quality
    ambiguities: int       # overlap positions where an IUPAC code was emitted
    warnings: list[str]

    @property
    def total_length(self) -> int:
        return len(self.consensus_seq)

    @property
    def overlap_length(self) -> int:
        return self.overlap.overlap_length


# ── Core functions ────────────────────────────────────────────────────────────

def reverse_complement_with_quality(
    seq: str, quality: np.ndarray
) -> tuple[str, np.ndarray]:
    """Return RC of *seq* and the quality array in reversed order."""
    return reverse_complement(seq), quality[::-1].copy()


def _extract_aligned_strings(
    seq1: str, seq2: str, coords: np.ndarray
) -> tuple[str, str]:
    """Build gapped alignment strings from Bio.Align coordinate array.

    coords: shape (2, n) from PairwiseAlignment.coordinates.
    coords[0] = positions in seq1 (target); coords[1] = positions in seq2 (query).
    """
    parts1: list[str] = []
    parts2: list[str] = []
    t = coords[0]
    q = coords[1]
    for i in range(len(t) - 1):
        dt = int(t[i + 1]) - int(t[i])
        dq = int(q[i + 1]) - int(q[i])
        if dt > 0 and dq > 0:
            # Diagonal block: match or mismatch (dt == dq always holds here)
            parts1.append(seq1[t[i] : t[i + 1]])
            parts2.append(seq2[q[i] : q[i + 1]])
        elif dt > 0:
            # Gap in query (seq2)
            parts1.append(seq1[t[i] : t[i + 1]])
            parts2.append("-" * dt)
        else:
            # Gap in target (seq1)
            parts1.append("-" * dq)
            parts2.append(seq2[q[i] : q[i + 1]])
    return "".join(parts1), "".join(parts2)


def find_overlap(
    fwd_seq: str,
    rc_seq: str,
    min_overlap: int = 20,
    min_identity: float = 0.90,
    match_score: float = 2.0,
    mismatch_score: float = -1.0,
    open_gap_score: float = -2.5,
    extend_gap_score: float = -0.5,
) -> OverlapResult:
    """Smith-Waterman local alignment to locate the overlap region.

    Returns an OverlapResult; ``has_overlap`` is False when the best alignment
    score falls below the threshold or the aligned region is shorter than
    *min_overlap*.
    """
    _no_overlap = OverlapResult(
        fwd_start=0, fwd_end=0, rc_start=0, rc_end=0,
        fwd_aligned="", rc_aligned="",
        overlap_length=0, identity=0.0, score=0.0, has_overlap=False,
    )

    if not fwd_seq or not rc_seq:
        return _no_overlap

    from Bio import Align

    aligner = Align.PairwiseAligner(mode="local")
    aligner.match_score = match_score
    aligner.mismatch_score = mismatch_score
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = extend_gap_score

    try:
        alignments = aligner.align(fwd_seq.upper(), rc_seq.upper())
        best = next(iter(alignments))
    except StopIteration:
        return _no_overlap

    score = float(best.score)
    min_score = min_overlap * match_score * 0.7  # require ≥70% match rate
    if score <= 0:
        return _no_overlap

    coords = best.coordinates  # shape (2, n)
    fwd_start = int(coords[0][0])
    fwd_end   = int(coords[0][-1])
    rc_start  = int(coords[1][0])
    rc_end    = int(coords[1][-1])

    fwd_aln, rc_aln = _extract_aligned_strings(fwd_seq, rc_seq, coords)
    overlap_len = len(fwd_aln)

    # Count identity over non-gap columns
    matched = sum(
        1 for f, r in zip(fwd_aln, rc_aln)
        if f != "-" and r != "-" and f.upper() == r.upper()
    )
    total_cols = sum(1 for f, r in zip(fwd_aln, rc_aln) if f != "-" and r != "-")
    identity = matched / total_cols if total_cols > 0 else 0.0

    actual_span = min(fwd_end - fwd_start, rc_end - rc_start)
    has_overlap = (
        score >= min_score
        and actual_span >= min_overlap
        and identity >= min_identity
    )

    return OverlapResult(
        fwd_start=fwd_start,
        fwd_end=fwd_end,
        rc_start=rc_start,
        rc_end=rc_end,
        fwd_aligned=fwd_aln,
        rc_aligned=rc_aln,
        overlap_length=overlap_len,
        identity=identity,
        score=score,
        has_overlap=has_overlap,
    )


def build_consensus(
    fwd_seq: str,
    fwd_quality: np.ndarray,
    rc_seq: str,
    rc_quality: np.ndarray,
    overlap: OverlapResult,
    low_quality_threshold: int = 20,
    quality_diff_threshold: int = 10,
) -> ConsensusResult:
    """Build a per-position consensus from forward and RC-reverse reads.

    Outside the overlap: take the covering read verbatim.
    Inside the overlap:
      - Agree → base, quality = min(60, q_fwd + q_rc).
      - Disagree, Δq ≥ *quality_diff_threshold* → higher-quality base.
      - Disagree, Δq < threshold or both low → IUPAC ambiguity code.
      - One side has gap (-) → use the other side's base.
    """
    positions: list[ConsensusPosition] = []
    disagreements = 0
    ambiguities = 0
    warnings: list[str] = []

    def _q(arr: np.ndarray, i: int) -> int:
        return int(arr[i]) if 0 <= i < len(arr) else 0

    # ── Fwd-only prefix ──────────────────────────────────────────────────────
    for i in range(overlap.fwd_start):
        positions.append(ConsensusPosition(
            base=fwd_seq[i].upper(), quality=_q(fwd_quality, i), source="fwd"
        ))

    # ── Overlap region ───────────────────────────────────────────────────────
    fi = overlap.fwd_start
    ri = overlap.rc_start
    for fc, rcc in zip(overlap.fwd_aligned, overlap.rc_aligned):
        qf = _q(fwd_quality, fi)
        qr = _q(rc_quality, ri)

        if fc == "-":
            # Insertion from rc
            positions.append(ConsensusPosition(
                base=rcc.upper(), quality=qr, source="rev"
            ))
            ri += 1
        elif rcc == "-":
            # Insertion from fwd
            positions.append(ConsensusPosition(
                base=fc.upper(), quality=qf, source="fwd"
            ))
            fi += 1
        else:
            fbu = fc.upper()
            rbu = rcc.upper()
            if fbu == rbu:
                positions.append(ConsensusPosition(
                    base=fbu, quality=min(60, qf + qr), source="agree"
                ))
            else:
                qdiff = abs(qf - qr)
                if qdiff >= quality_diff_threshold:
                    winner_base = fbu if qf >= qr else rbu
                    winner_q    = max(qf, qr)
                    positions.append(ConsensusPosition(
                        base=winner_base, quality=winner_q, source="disagree"
                    ))
                    disagreements += 1
                else:
                    positions.append(ConsensusPosition(
                        base=_iupac(fbu, rbu),
                        quality=min(qf, qr),
                        source="ambiguous",
                    ))
                    ambiguities += 1
            fi += 1
            ri += 1

    # ── RC-only suffix ────────────────────────────────────────────────────────
    for i in range(overlap.rc_end, len(rc_seq)):
        positions.append(ConsensusPosition(
            base=rc_seq[i].upper(), quality=_q(rc_quality, i), source="rev"
        ))

    # ── Sanity warnings ───────────────────────────────────────────────────────
    if overlap.fwd_end < len(fwd_seq):
        n = len(fwd_seq) - overlap.fwd_end
        warnings.append(
            f"{n} base{'s' if n != 1 else ''} at the 3' end of the forward read "
            "lie outside the overlap and are not included in the consensus."
        )
    if overlap.rc_start > 0:
        n = overlap.rc_start
        warnings.append(
            f"{n} base{'s' if n != 1 else ''} at the 5' end of the RC reverse read "
            "lie outside the overlap and are not included in the consensus."
        )

    return ConsensusResult(
        consensus_seq="".join(p.base for p in positions),
        consensus_quality=[p.quality for p in positions],
        positions=positions,
        overlap=overlap,
        fwd_only_length=overlap.fwd_start,
        rev_only_length=max(0, len(rc_seq) - overlap.rc_end),
        disagreements=disagreements,
        ambiguities=ambiguities,
        warnings=warnings,
    )
