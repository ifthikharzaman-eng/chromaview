"""Tests for quality analysis module."""
import numpy as np
import pytest

from chromaview.core.models import Chromatogram
from chromaview.analysis.quality import trim_low_quality, quality_summary
from chromaview.tests.test_models import make_test_chrom


class TestTrimming:
    def test_basic_trim(self):
        chrom = make_test_chrom(n_bases=100)
        # Set low quality at edges, high in middle
        chrom.quality_scores[:15] = 5
        chrom.quality_scores[-15:] = 5
        chrom.quality_scores[15:85] = 40

        start, end = trim_low_quality(chrom, threshold=20)
        assert start >= 10
        assert end <= 90
        assert end > start

    def test_all_high_quality(self):
        chrom = make_test_chrom(n_bases=50)
        chrom.quality_scores[:] = 40
        start, end = trim_low_quality(chrom, threshold=20)
        assert start == 0
        assert end == 50

    def test_all_low_quality(self):
        chrom = make_test_chrom(n_bases=50)
        chrom.quality_scores[:] = 5
        start, end = trim_low_quality(chrom, threshold=20)
        # Should fall back to full range when nothing passes
        assert end - start > 0


class TestQualitySummary:
    def test_basic_stats(self):
        chrom = make_test_chrom(n_bases=100)
        chrom.quality_scores[:50] = 35
        chrom.quality_scores[50:] = 15
        qs = quality_summary(chrom)

        assert qs.total_bases == 100
        assert qs.q20_count > 0
        assert qs.q30_count > 0
        assert 0 <= qs.q20_fraction <= 1
        assert 0 <= qs.q30_fraction <= 1
        assert qs.min_quality >= 0
        assert qs.max_quality <= 60

    def test_empty_chrom(self):
        chrom = Chromatogram()
        qs = quality_summary(chrom)
        assert qs.total_bases == 0
        assert qs.mean_quality == 0.0
