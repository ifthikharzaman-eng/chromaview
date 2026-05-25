"""Tests for the Chromatogram data model."""
import numpy as np
import pytest

from chromaview.core.models import Chromatogram, ChromatogramMetadata


def make_test_chrom(n_bases=100, n_samples=5000):
    """Create a synthetic chromatogram for testing."""
    rng = np.random.default_rng(42)
    bases = list(rng.choice(list("ACGT"), size=n_bases))
    peak_locations = np.linspace(100, n_samples - 100, n_bases).astype(int)
    quality = rng.integers(5, 55, size=n_bases)

    traces = {}
    for base in "ACGT":
        trace = rng.normal(50, 20, size=n_samples).clip(0)
        # Add peaks at called positions
        for i, (b, loc) in enumerate(zip(bases, peak_locations)):
            if b == base:
                trace[max(0, loc - 5) : loc + 5] += 500
            else:
                trace[max(0, loc - 5) : loc + 5] += rng.uniform(10, 50)
        traces[base] = trace

    return Chromatogram(
        traces=traces,
        raw_traces={},
        basecalls=bases,
        peak_locations=peak_locations,
        quality_scores=quality,
        metadata=ChromatogramMetadata(sample_name="test_sample", file_format="ab1"),
    )


class TestChromatogram:
    def test_basic_properties(self):
        chrom = make_test_chrom(n_bases=50)
        assert chrom.num_bases == 50
        assert chrom.trace_length == 5000
        assert len(chrom.sequence) == 50

    def test_trimming(self):
        chrom = make_test_chrom(n_bases=100)
        chrom.trim_start = 10
        chrom.trim_end = 90
        assert len(chrom.trimmed_basecalls) == 80
        assert len(chrom.trimmed_quality) == 80
        assert len(chrom.trimmed_sequence) == 80

    def test_edit_base(self):
        chrom = make_test_chrom()
        original = chrom.basecalls[5]
        new_base = "T" if original != "T" else "A"
        chrom.edit_base(5, new_base)
        assert chrom.basecalls[5] == new_base
        assert chrom.is_edited
        assert len(chrom.get_edits()) == 1

    def test_undo_redo(self):
        chrom = make_test_chrom()
        original = chrom.basecalls[5]
        chrom.edit_base(5, "N")
        assert chrom.basecalls[5] == "N"

        chrom.undo()
        assert chrom.basecalls[5] == original

        chrom.redo()
        assert chrom.basecalls[5] == "N"

    def test_reset_edits(self):
        chrom = make_test_chrom()
        original_seq = chrom.sequence
        chrom.edit_base(0, "N")
        chrom.edit_base(1, "N")
        chrom.reset_edits()
        assert chrom.sequence == original_seq
        assert not chrom.is_edited

    def test_edit_invalid_base(self):
        chrom = make_test_chrom()
        with pytest.raises(ValueError):
            chrom.edit_base(0, "X")

    def test_edit_out_of_range(self):
        chrom = make_test_chrom(n_bases=10)
        with pytest.raises(IndexError):
            chrom.edit_base(100, "A")

    def test_mean_quality(self):
        chrom = make_test_chrom()
        assert chrom.mean_quality > 0

    def test_clone(self):
        chrom = make_test_chrom()
        clone = chrom.clone()
        clone.edit_base(0, "N")
        assert chrom.basecalls[0] != "N"  # original unchanged
