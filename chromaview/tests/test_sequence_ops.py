"""Tests for sequence operations module."""
import pytest

from chromaview.core.sequence_ops import (
    reverse_complement,
    align_to_reference,
    search_sequence,
)


class TestReverseComplement:
    def test_basic(self):
        assert reverse_complement("ATGC") == "GCAT"

    def test_palindrome(self):
        assert reverse_complement("AATT") == "AATT"

    def test_ambiguity_codes(self):
        assert reverse_complement("RYSWKM") == "KMWSRY"

    def test_empty(self):
        assert reverse_complement("") == ""

    def test_single_base(self):
        assert reverse_complement("A") == "T"
        assert reverse_complement("G") == "C"

    def test_gap(self):
        assert reverse_complement("AT-GC") == "GC-AT"


class TestAlignment:
    def test_identical(self):
        result = align_to_reference("ATGCATGC", "ATGCATGC")
        assert result.identity == 1.0
        assert len(result.mismatch_positions) == 0

    def test_single_snp(self):
        result = align_to_reference("ATGCATGC", "ATGAATGC")
        assert result.identity < 1.0
        assert len(result.mismatch_positions) >= 1

    def test_empty(self):
        result = align_to_reference("", "ATGC")
        assert result.identity == 0.0

    def test_different_lengths(self):
        result = align_to_reference("ATGC", "ATGCATGC")
        assert result.identity > 0  # aligned portion should match
        assert len(result.query_aligned) == len(result.ref_aligned)


class TestSearch:
    def test_motif_found(self):
        hits = search_sequence("ATGCGATCGATCG", "GATC")
        assert len(hits) >= 2

    def test_motif_not_found(self):
        hits = search_sequence("AAAA", "CCCC")
        assert len(hits) == 0

    def test_position_search(self):
        hits = search_sequence("ATGCATGC", "3")
        assert len(hits) == 1
        assert hits[0].start == 2  # 0-based

    def test_case_insensitive(self):
        hits = search_sequence("ATGC", "atgc")
        assert len(hits) == 1

    def test_position_out_of_range(self):
        hits = search_sequence("ATG", "100")
        assert len(hits) == 0
