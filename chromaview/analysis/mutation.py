"""Basic mutation calling and ambiguous base detection.

Detects positions where multiple channels show significant signal,
suggesting heterozygous mutations or mixed bases. Also provides
simple SNP calling when compared against a reference.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core.models import Chromatogram


@dataclass
class AmbiguousBase:
    """A position where the chromatogram shows potential mixed signal."""

    position: int           # 0-based index in basecalls
    called_base: str
    suggested_iupac: str    # IUPAC ambiguity code
    primary_height: float
    secondary_height: float
    secondary_base: str
    ratio: float            # secondary / primary


# IUPAC codes for two-base ambiguities
_IUPAC_PAIR = {
    frozenset({"A", "G"}): "R",
    frozenset({"C", "T"}): "Y",
    frozenset({"G", "C"}): "S",
    frozenset({"A", "T"}): "W",
    frozenset({"G", "T"}): "K",
    frozenset({"A", "C"}): "M",
}


def detect_ambiguous_bases(
    chrom: Chromatogram,
    min_ratio: float = 0.30,
    min_secondary_height: float = 100.0,
) -> list[AmbiguousBase]:
    """Find positions with significant secondary peaks.

    A secondary peak is considered significant if:
        - Its height is >= min_ratio × the primary peak height
        - Its absolute height is >= min_secondary_height

    Args:
        chrom: Chromatogram to analyze.
        min_ratio: Minimum secondary/primary ratio to flag (default 0.30 = 30%).
        min_secondary_height: Absolute minimum for the secondary peak.

    Returns:
        List of AmbiguousBase instances for flagged positions.
    """
    all_bases = ["A", "T", "G", "C"]
    results = []

    for i, (base, loc) in enumerate(zip(chrom.basecalls, chrom.peak_locations)):
        loc = int(loc)
        base = base.upper()

        # Get heights for all four channels at this peak
        heights = {}
        for b in all_bases:
            trace = chrom.traces.get(b)
            if trace is not None and loc < len(trace):
                heights[b] = float(trace[loc])
            else:
                heights[b] = 0.0

        primary = heights.get(base, 0.0)
        if primary <= 0:
            continue

        # Find strongest secondary channel
        secondary_bases = [(b, h) for b, h in heights.items() if b != base]
        secondary_bases.sort(key=lambda x: x[1], reverse=True)
        sec_base, sec_height = secondary_bases[0]

        ratio = sec_height / primary if primary > 0 else 0.0

        if ratio >= min_ratio and sec_height >= min_secondary_height:
            pair = frozenset({base, sec_base})
            iupac = _IUPAC_PAIR.get(pair, "N")
            results.append(AmbiguousBase(
                position=i,
                called_base=base,
                suggested_iupac=iupac,
                primary_height=primary,
                secondary_height=sec_height,
                secondary_base=sec_base,
                ratio=ratio,
            ))

    return results


@dataclass
class MutationCall:
    """A potential mutation found by comparing to a reference."""

    query_position: int    # 0-based in query
    ref_position: int      # 0-based in reference
    query_base: str
    ref_base: str
    quality: int           # Phred score at this position
    mutation_type: str     # "SNP", "insertion", "deletion"


def call_mutations(
    query: str,
    reference: str,
    quality_scores: np.ndarray | None = None,
    alignment_result=None,
) -> list[MutationCall]:
    """Compare query against reference and report mutations.

    If an AlignmentResult is provided, uses the pre-computed alignment.
    Otherwise does a simple positional comparison (no gaps).
    """
    mutations = []

    if alignment_result is not None:
        q_aln = alignment_result.query_aligned
        r_aln = alignment_result.ref_aligned
        q_pos = -1
        r_pos = -1

        for aln_pos, (qb, rb) in enumerate(zip(q_aln, r_aln)):
            if qb != "-":
                q_pos += 1
            if rb != "-":
                r_pos += 1

            if qb == "-" and rb != "-":
                mutations.append(MutationCall(
                    query_position=q_pos,
                    ref_position=r_pos,
                    query_base="-",
                    ref_base=rb,
                    quality=0,
                    mutation_type="deletion",
                ))
            elif qb != "-" and rb == "-":
                q_qual = int(quality_scores[q_pos]) if quality_scores is not None and q_pos < len(quality_scores) else 0
                mutations.append(MutationCall(
                    query_position=q_pos,
                    ref_position=r_pos,
                    query_base=qb,
                    ref_base="-",
                    quality=q_qual,
                    mutation_type="insertion",
                ))
            elif qb.upper() != rb.upper():
                q_qual = int(quality_scores[q_pos]) if quality_scores is not None and q_pos < len(quality_scores) else 0
                mutations.append(MutationCall(
                    query_position=q_pos,
                    ref_position=r_pos,
                    query_base=qb,
                    ref_base=rb,
                    quality=q_qual,
                    mutation_type="SNP",
                ))
    else:
        # Simple positional comparison (no alignment)
        for i, (qb, rb) in enumerate(zip(query.upper(), reference.upper())):
            if qb != rb:
                q_qual = int(quality_scores[i]) if quality_scores is not None and i < len(quality_scores) else 0
                mutations.append(MutationCall(
                    query_position=i,
                    ref_position=i,
                    query_base=qb,
                    ref_base=rb,
                    quality=q_qual,
                    mutation_type="SNP",
                ))

    return mutations
