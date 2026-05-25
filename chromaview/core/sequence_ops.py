"""Sequence operations for chromatogram analysis.

Provides reverse complement, simple pairwise alignment against a reference,
and motif/position search within a called sequence.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


# IUPAC complement table (includes ambiguity codes)
_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn-",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn-",
)


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence string."""
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass
class AlignmentResult:
    """Result of aligning a query sequence to a reference."""

    query_aligned: str
    ref_aligned: str
    matches: list[bool]       # True where bases match
    mismatch_positions: list[int]  # 0-based positions in query
    identity: float           # fraction of matching bases
    score: int
    offset: int               # position in reference where query begins


def align_to_reference(
    query: str,
    reference: str,
    match_score: int = 2,
    mismatch_penalty: int = -1,
    gap_penalty: int = -2,
) -> AlignmentResult:
    """Simple Needleman-Wunsch global alignment of query against reference.

    For production use consider using parasail or Bio.Align, but this
    standalone implementation avoids heavy dependencies and is sufficient
    for single-read comparisons (< 1500 bp).
    """
    m, n = len(query), len(reference)
    if m == 0 or n == 0:
        return AlignmentResult("", "", [], [], 0.0, 0, 0)

    # ── Score matrix ─────────────────────────────────────────────
    dp = np.zeros((m + 1, n + 1), dtype=int)
    for i in range(1, m + 1):
        dp[i][0] = dp[i - 1][0] + gap_penalty
    for j in range(1, n + 1):
        dp[0][j] = dp[0][j - 1] + gap_penalty

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if query[i - 1].upper() == reference[j - 1].upper():
                diag = dp[i - 1][j - 1] + match_score
            else:
                diag = dp[i - 1][j - 1] + mismatch_penalty
            up = dp[i - 1][j] + gap_penalty
            left = dp[i][j - 1] + gap_penalty
            dp[i][j] = max(diag, up, left)

    # ── Traceback ────────────────────────────────────────────────
    q_aln, r_aln = [], []
    i, j = m, n
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            if query[i - 1].upper() == reference[j - 1].upper():
                s = match_score
            else:
                s = mismatch_penalty
            if dp[i][j] == dp[i - 1][j - 1] + s:
                q_aln.append(query[i - 1])
                r_aln.append(reference[j - 1])
                i -= 1
                j -= 1
                continue
        if i > 0 and dp[i][j] == dp[i - 1][j] + gap_penalty:
            q_aln.append(query[i - 1])
            r_aln.append("-")
            i -= 1
        else:
            q_aln.append("-")
            r_aln.append(reference[j - 1])
            j -= 1

    q_aln.reverse()
    r_aln.reverse()
    query_aligned = "".join(q_aln)
    ref_aligned = "".join(r_aln)

    # ── Compute matches ──────────────────────────────────────────
    matches = [
        q.upper() == r.upper() and q != "-" and r != "-"
        for q, r in zip(query_aligned, ref_aligned)
    ]
    mismatch_positions = [
        i for i, (q, r) in enumerate(zip(query_aligned, ref_aligned))
        if q.upper() != r.upper() and q != "-" and r != "-"
    ]
    total_aligned = sum(1 for q, r in zip(query_aligned, ref_aligned) if q != "-" and r != "-")
    identity = sum(matches) / total_aligned if total_aligned > 0 else 0.0

    return AlignmentResult(
        query_aligned=query_aligned,
        ref_aligned=ref_aligned,
        matches=matches,
        mismatch_positions=mismatch_positions,
        identity=identity,
        score=int(dp[m][n]),
        offset=0,
    )


@dataclass
class SearchHit:
    """A single search result within a sequence."""
    start: int   # 0-based position in basecalls
    end: int     # exclusive end
    motif: str


def search_sequence(sequence: str, query: str) -> list[SearchHit]:
    """Find all occurrences of *query* in *sequence* (case-insensitive).

    Also supports position lookup: if query is a number, returns
    a single hit at that 1-based position.
    """
    hits = []
    # Position-based search
    if query.isdigit():
        pos = int(query) - 1  # convert 1-based to 0-based
        if 0 <= pos < len(sequence):
            hits.append(SearchHit(start=pos, end=pos + 1, motif=sequence[pos]))
        return hits

    # Motif search
    seq_upper = sequence.upper()
    q_upper = query.upper()
    start = 0
    while True:
        idx = seq_upper.find(q_upper, start)
        if idx == -1:
            break
        hits.append(SearchHit(start=idx, end=idx + len(q_upper), motif=query))
        start = idx + 1
    return hits
