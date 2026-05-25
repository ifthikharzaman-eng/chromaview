"""Tests for paired-read consensus assembly.

All tests use synthetic sequences — no file I/O, no network.
Covers:
  - reverse_complement_with_quality
  - find_overlap: perfect overlap, no overlap, short overlap, mismatch
  - build_consensus: perfect overlap, IUPAC on equal-quality mismatch,
    quality winner on high-Δq mismatch, fwd-tail warning, rc-prefix warning

Sequence design rationale
--------------------------
  _FWD_SEQ = "A"*20  +  GC-repeat-25  (45 bp)
  _RC_SEQ  =             GC-repeat-25  +  "T"*20  (45 bp)

The three regions have distinct base compositions: A-only prefix, GC-only
overlap, T-only suffix.  Because A/T do not match G/C, Smith-Waterman cannot
extend the local alignment beyond the GC block, giving unambiguous boundaries:
  fwd_start=20, fwd_end=45, rc_start=0, rc_end=25.
"""
from __future__ import annotations

import numpy as np
import pytest

from chromaview.analysis.consensus import (
    ConsensusResult,
    OverlapResult,
    _iupac,
    build_consensus,
    find_overlap,
    reverse_complement_with_quality,
)

# ── Shared test data ──────────────────────────────────────────────────────────

# 25-bp GC-only overlap; A/T regions cannot extend the SW alignment into it.
_OVERLAP = "GC" * 12 + "G"          # "GCGCGCGCGCGCGCGCGCGCGCGCG", 25 bp
assert len(_OVERLAP) == 25

_FWD_SEQ = "A" * 20 + _OVERLAP      # 45 bp
_RC_SEQ  = _OVERLAP + "T" * 20      # 45 bp

assert _FWD_SEQ[20:] == _OVERLAP    # last 25 of fwd == overlap
assert _RC_SEQ[:25]  == _OVERLAP    # first 25 of rc == overlap

_Q40 = np.full(45, 40, dtype=int)   # uniform Phred-40 quality

# One-mismatch variant: change _OVERLAP[12] ('G') → 'A' in the rc version
# IUPAC for G+A = R;  _OVERLAP[12] is index 12 (even → 'G')
_OVERLAP_MUT = _OVERLAP[:12] + "A" + _OVERLAP[13:]   # 25 bp
_RC_SEQ_MUT  = _OVERLAP_MUT + "T" * 20               # 45 bp


# ── _iupac helper ─────────────────────────────────────────────────────────────

class TestIUPAC:
    def test_ac_is_M(self):
        assert _iupac("A", "C") == "M"
        assert _iupac("C", "A") == "M"

    def test_gt_is_K(self):
        assert _iupac("G", "T") == "K"

    def test_ag_is_R(self):
        assert _iupac("A", "G") == "R"
        assert _iupac("G", "A") == "R"

    def test_case_insensitive(self):
        assert _iupac("a", "g") == "R"

    def test_unknown_pair_returns_N(self):
        assert _iupac("A", "N") == "N"


# ── reverse_complement_with_quality ──────────────────────────────────────────

class TestReverseComplementWithQuality:
    def test_sequence_is_rc(self):
        seq, _qual = reverse_complement_with_quality("ATCG", np.array([10, 20, 30, 40]))
        assert seq == "CGAT"

    def test_quality_is_reversed(self):
        _seq, qual = reverse_complement_with_quality("ATCG", np.array([10, 20, 30, 40]))
        assert list(qual) == [40, 30, 20, 10]

    def test_does_not_modify_original(self):
        orig = np.array([1, 2, 3, 4])
        _, rc_q = reverse_complement_with_quality("ACGT", orig)
        rc_q[0] = 99
        assert orig[0] == 1


# ── find_overlap ─────────────────────────────────────────────────────────────

class TestFindOverlap:
    def test_perfect_overlap_detected(self):
        result = find_overlap(_FWD_SEQ, _RC_SEQ)
        assert result.has_overlap

    def test_perfect_overlap_fwd_boundaries(self):
        result = find_overlap(_FWD_SEQ, _RC_SEQ)
        # GC block starts at index 20 in fwd and ends at 45 (exclusive)
        assert result.fwd_start == 20
        assert result.fwd_end == 45

    def test_perfect_overlap_rc_boundaries(self):
        result = find_overlap(_FWD_SEQ, _RC_SEQ)
        assert result.rc_start == 0
        assert result.rc_end == 25

    def test_perfect_overlap_identity_is_1(self):
        result = find_overlap(_FWD_SEQ, _RC_SEQ)
        assert result.identity == pytest.approx(1.0)

    def test_perfect_overlap_length_is_25(self):
        result = find_overlap(_FWD_SEQ, _RC_SEQ)
        assert result.overlap_length == 25

    def test_perfect_score_positive(self):
        result = find_overlap(_FWD_SEQ, _RC_SEQ)
        assert result.score > 0

    def test_no_overlap_all_mismatch(self):
        # A vs C: every aligned pair is a mismatch → SW score 0 → no overlap
        result = find_overlap("A" * 30, "C" * 30)
        assert not result.has_overlap

    def test_empty_fwd_returns_false(self):
        assert not find_overlap("", "ACGT").has_overlap

    def test_empty_rc_returns_false(self):
        assert not find_overlap("ACGT", "").has_overlap

    def test_overlap_below_min_returns_false(self):
        # 10-bp overlap with min_overlap=20 → rejected
        fwd = "A" * 20 + "GC" * 5           # overlap = "GCGCGCGCGC" (10 bp)
        rc  = "GC" * 5 + "T" * 20
        result = find_overlap(fwd, rc, min_overlap=20)
        assert not result.has_overlap

    def test_mismatch_identity_below_1(self):
        # One mismatch in 25-bp overlap → identity = 24/25
        result = find_overlap(_FWD_SEQ, _RC_SEQ_MUT)
        assert result.has_overlap
        assert result.identity == pytest.approx(24 / 25)

    def test_mismatch_fwd_boundaries_unchanged(self):
        result = find_overlap(_FWD_SEQ, _RC_SEQ_MUT)
        assert result.fwd_start == 20
        assert result.fwd_end == 45

    def test_identity_below_threshold_rejected(self):
        # min_identity=0.99 rejects the one-mismatch overlap (identity=24/25=0.96)
        result = find_overlap(_FWD_SEQ, _RC_SEQ_MUT, min_identity=0.99)
        assert not result.has_overlap

    def test_identity_threshold_identity_still_reported(self):
        # identity is measured and returned even when has_overlap=False
        result = find_overlap(_FWD_SEQ, _RC_SEQ_MUT, min_identity=0.99)
        assert result.identity == pytest.approx(24 / 25)


# ── build_consensus ───────────────────────────────────────────────────────────

class TestBuildConsensus:
    # ── perfect overlap ──────────────────────────────────────────────────────

    def test_perfect_consensus_length(self):
        ov = find_overlap(_FWD_SEQ, _RC_SEQ)
        result = build_consensus(_FWD_SEQ, _Q40, _RC_SEQ, _Q40, ov)
        assert result.total_length == 65   # 20 + 25 + 20

    def test_perfect_consensus_sequence(self):
        ov = find_overlap(_FWD_SEQ, _RC_SEQ)
        result = build_consensus(_FWD_SEQ, _Q40, _RC_SEQ, _Q40, ov)
        expected = "A" * 20 + _OVERLAP + "T" * 20
        assert result.consensus_seq == expected

    def test_perfect_no_disagreements_or_ambiguities(self):
        ov = find_overlap(_FWD_SEQ, _RC_SEQ)
        result = build_consensus(_FWD_SEQ, _Q40, _RC_SEQ, _Q40, ov)
        assert result.disagreements == 0
        assert result.ambiguities == 0

    def test_agree_positions_quality_combined(self):
        ov = find_overlap(_FWD_SEQ, _RC_SEQ)
        result = build_consensus(_FWD_SEQ, _Q40, _RC_SEQ, _Q40, ov)
        # quality = min(60, 40+40) = 60 in the overlap
        agree_quals = [p.quality for p in result.positions if p.source == "agree"]
        assert all(q == 60 for q in agree_quals)

    def test_fwd_only_prefix_length(self):
        ov = find_overlap(_FWD_SEQ, _RC_SEQ)
        result = build_consensus(_FWD_SEQ, _Q40, _RC_SEQ, _Q40, ov)
        assert result.fwd_only_length == 20

    def test_rev_only_suffix_length(self):
        ov = find_overlap(_FWD_SEQ, _RC_SEQ)
        result = build_consensus(_FWD_SEQ, _Q40, _RC_SEQ, _Q40, ov)
        assert result.rev_only_length == 20

    def test_no_warnings_on_clean_overlap(self):
        ov = find_overlap(_FWD_SEQ, _RC_SEQ)
        result = build_consensus(_FWD_SEQ, _Q40, _RC_SEQ, _Q40, ov)
        assert result.warnings == []

    # ── one mismatch → IUPAC (equal quality) ────────────────────────────────

    def test_mismatch_equal_quality_emits_one_iupac(self):
        ov = find_overlap(_FWD_SEQ, _RC_SEQ_MUT)
        result = build_consensus(_FWD_SEQ, _Q40, _RC_SEQ_MUT, _Q40, ov)
        assert result.ambiguities == 1
        assert result.disagreements == 0

    def test_mismatch_equal_quality_iupac_code_is_R(self):
        # fwd has 'G' at overlap[12], rc_mut has 'A' → G+A = R
        ov = find_overlap(_FWD_SEQ, _RC_SEQ_MUT)
        result = build_consensus(_FWD_SEQ, _Q40, _RC_SEQ_MUT, _Q40, ov)
        ambig = [p for p in result.positions if p.source == "ambiguous"]
        assert len(ambig) == 1
        assert ambig[0].base == "R"

    # ── one mismatch → higher-quality winner ────────────────────────────────

    def test_mismatch_quality_winner_is_fwd(self):
        # fwd q=40, rc q=10 at mismatch → Δq=30 ≥ threshold=10 → fwd wins
        q_rc_low = np.full(45, 10, dtype=int)
        ov = find_overlap(_FWD_SEQ, _RC_SEQ_MUT)
        result = build_consensus(_FWD_SEQ, _Q40, _RC_SEQ_MUT, q_rc_low, ov)
        assert result.disagreements == 1
        assert result.ambiguities == 0

    def test_mismatch_quality_winner_base_is_G(self):
        q_rc_low = np.full(45, 10, dtype=int)
        ov = find_overlap(_FWD_SEQ, _RC_SEQ_MUT)
        result = build_consensus(_FWD_SEQ, _Q40, _RC_SEQ_MUT, q_rc_low, ov)
        disagree = [p for p in result.positions if p.source == "disagree"]
        assert len(disagree) == 1
        assert disagree[0].base == "G"   # fwd (q=40) beats rc (q=10)

    # ── warnings ─────────────────────────────────────────────────────────────

    def test_fwd_tail_warning(self):
        # fwd_end < len(fwd) → warning about unused 3′ tail
        fake_ov = OverlapResult(
            fwd_start=5, fwd_end=20,          # 25 fwd chars after the overlap
            rc_start=0,  rc_end=15,
            fwd_aligned=_FWD_SEQ[5:20],       # 15 A's
            rc_aligned=_RC_SEQ[0:15],         # 15 GC chars
            overlap_length=15,
            identity=0.0, score=30.0, has_overlap=True,
        )
        q = np.full(45, 30, dtype=int)
        result = build_consensus(_FWD_SEQ, q, _RC_SEQ, q, fake_ov)
        assert any("3'" in w for w in result.warnings)

    def test_rc_prefix_warning(self):
        # rc_start > 0 → warning about unused 5′ rc prefix
        fake_ov = OverlapResult(
            fwd_start=0, fwd_end=20,
            rc_start=5,  rc_end=25,           # 5 rc chars before the overlap
            fwd_aligned=_FWD_SEQ[0:20],       # 20 A's
            rc_aligned=_RC_SEQ[5:25],         # 20 GC chars
            overlap_length=20,
            identity=0.0, score=40.0, has_overlap=True,
        )
        q = np.full(45, 30, dtype=int)
        result = build_consensus(_FWD_SEQ, q, _RC_SEQ, q, fake_ov)
        assert any("5'" in w for w in result.warnings)

    # ── ConsensusResult properties ───────────────────────────────────────────

    def test_total_length_equals_seq_length(self):
        ov = find_overlap(_FWD_SEQ, _RC_SEQ)
        result = build_consensus(_FWD_SEQ, _Q40, _RC_SEQ, _Q40, ov)
        assert result.total_length == len(result.consensus_seq)

    def test_overlap_length_property_delegates(self):
        ov = find_overlap(_FWD_SEQ, _RC_SEQ)
        result = build_consensus(_FWD_SEQ, _Q40, _RC_SEQ, _Q40, ov)
        assert result.overlap_length == ov.overlap_length
