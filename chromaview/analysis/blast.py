"""BLAST search integration for ChromaView.

Pure data layer: BLASTHit dataclass + parse_blast_record().
Network I/O lives in gui/blast_dialog.py (BLASTWorker) so this
module stays importable without PyQt6.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

MIN_QUERY_LENGTH = 30  # warn if trimmed sequence is shorter

DATABASES = [
    ("core_nt", "core_nt — Core nucleotide (default)"),
    ("nt", "nt — All nucleotide"),
    ("refseq_rna", "refseq_rna — RefSeq RNA"),
]

_ORG_RE = re.compile(r"\[([^\[\]]+)\]\s*$")


@dataclass
class BLASTHit:
    """One BLAST hit (best HSP only)."""
    accession: str
    description: str
    organism: str
    query_length: int
    hsp_query_start: int   # 0-based, inclusive
    hsp_query_end: int     # 0-based, exclusive
    hsp_align_length: int  # alignment columns (including gaps)
    identity_count: int    # number of identical positions
    evalue: float
    bit_score: float
    score: float
    alignment_text: str

    @property
    def query_coverage(self) -> float:
        if self.query_length <= 0:
            return 0.0
        return (self.hsp_query_end - self.hsp_query_start) / self.query_length

    @property
    def pct_identity(self) -> float:
        if self.hsp_align_length <= 0:
            return 0.0
        return self.identity_count / self.hsp_align_length


def _extract_organism(description: str) -> str:
    m = _ORG_RE.search(description.strip())
    return m.group(1) if m else ""


def parse_blast_record(record) -> list[BLASTHit]:
    """Convert a Bio.Blast.Record to a list of BLASTHit objects.

    Uses the best (first) HSP per hit. Compatible with Bio.Blast 1.83+.
    """
    hits: list[BLASTHit] = []

    query_length = 0
    if hasattr(record, "query") and record.query is not None:
        try:
            query_length = len(record.query.seq)
        except Exception:
            pass

    for hit in record:
        if not hit:
            continue
        hsp = hit[0]

        accession = (hit.target.name or hit.target.id or "").strip()
        description = (hit.target.description or "").strip()
        organism = _extract_organism(description)

        # coordinates shape (2, n): row 0 = target, row 1 = query (0-based)
        coords = hsp.coordinates
        q_start = int(min(coords[1]))
        q_end = int(max(coords[1]))

        # hsp.shape == (nrows=2, ncols=alignment_length_including_gaps)
        align_length = int(hsp.shape[1])
        identity = int(hsp.annotations.get("identity", 0))
        evalue = float(hsp.annotations.get("evalue", float("inf")))
        bit_score = float(hsp.annotations.get("bit score", 0.0))
        score = float(hsp.score) if hsp.score is not None else 0.0
        alignment_text = str(hsp)

        hits.append(BLASTHit(
            accession=accession,
            description=description,
            organism=organism,
            query_length=query_length,
            hsp_query_start=q_start,
            hsp_query_end=q_end,
            hsp_align_length=align_length,
            identity_count=identity,
            evalue=evalue,
            bit_score=bit_score,
            score=score,
            alignment_text=alignment_text,
        ))

    return hits
