"""Tests for BLAST integration.

All network calls are mocked — tests run fully offline and fast.
Covers: XML parsing, hit extraction, organism extraction, query coverage,
percent identity, short-sequence warning threshold, and error handling paths.

Key design note: NCBIWWW.qblast() returns io.StringIO (text mode) because it
calls handle.read().decode() then wraps in StringIO.  Bio.Blast.read() requires
binary mode.  BLASTWorker bridges this by re-encoding to BytesIO.  The tests
below use StringIO as the mock return value to faithfully represent real behaviour.
"""
from __future__ import annotations

import io
from unittest.mock import patch, MagicMock

import pytest

from chromaview.analysis.blast import (
    BLASTHit,
    MIN_QUERY_LENGTH,
    DATABASES,
    parse_blast_record,
    _extract_organism,
)


# ── Minimal BLAST XML fixtures ────────────────────────────────────────────────

def _blast_xml(
    query_len: int = 100,
    hit_id: str = "ref|NM_000001.1|",
    hit_accession: str = "NM_000001",
    hit_def: str = "Homo sapiens test gene mRNA [Homo sapiens]",
    hit_len: int = 500,
    bit_score: float = 180.5,
    score: int = 97,
    evalue: str = "1.4e-45",
    q_from: int = 1,
    q_to: int = 97,
    h_from: int = 101,
    h_to: int = 197,
    identity: int = 95,
    positive: int = 95,
    gaps: int = 2,
    align_len: int = 97,
    n_hits: int = 1,
) -> bytes:
    seq = ("ATGC" * 25)[:align_len]
    mid = "|" * (align_len - 2) + "  "  # two mismatches at the end
    hits_xml = ""
    for i in range(n_hits):
        accession = f"{hit_accession}.{i}" if i > 0 else hit_accession
        hits_xml += f"""<Hit>
<Hit_num>{i + 1}</Hit_num>
<Hit_id>{hit_id}</Hit_id>
<Hit_def>{hit_def}</Hit_def>
<Hit_accession>{accession}</Hit_accession>
<Hit_len>{hit_len}</Hit_len>
<Hit_hsps><Hsp>
<Hsp_num>1</Hsp_num>
<Hsp_bit-score>{bit_score}</Hsp_bit-score>
<Hsp_score>{score}</Hsp_score>
<Hsp_evalue>{evalue}</Hsp_evalue>
<Hsp_query-from>{q_from}</Hsp_query-from>
<Hsp_query-to>{q_to}</Hsp_query-to>
<Hsp_hit-from>{h_from}</Hsp_hit-from>
<Hsp_hit-to>{h_to}</Hsp_hit-to>
<Hsp_query-frame>1</Hsp_query-frame>
<Hsp_hit-frame>1</Hsp_hit-frame>
<Hsp_identity>{identity}</Hsp_identity>
<Hsp_positive>{positive}</Hsp_positive>
<Hsp_gaps>{gaps}</Hsp_gaps>
<Hsp_align-len>{align_len}</Hsp_align-len>
<Hsp_qseq>{seq}</Hsp_qseq>
<Hsp_hseq>{seq}</Hsp_hseq>
<Hsp_midline>{mid}</Hsp_midline>
</Hsp></Hit_hsps>
</Hit>"""

    return f"""<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "NCBI_BlastOutput.dtd">
<BlastOutput>
<BlastOutput_program>blastn</BlastOutput_program>
<BlastOutput_version>BLASTN 2.14.0+</BlastOutput_version>
<BlastOutput_reference>Altschul SF et al.</BlastOutput_reference>
<BlastOutput_db>core_nt</BlastOutput_db>
<BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
<BlastOutput_query-def>TestQuery</BlastOutput_query-def>
<BlastOutput_query-len>{query_len}</BlastOutput_query-len>
<BlastOutput_param><Parameters>
<Parameters_expect>10</Parameters_expect>
<Parameters_sc-match>1</Parameters_sc-match>
<Parameters_sc-mismatch>-2</Parameters_sc-mismatch>
<Parameters_gap-open>0</Parameters_gap-open>
<Parameters_gap-extend>2</Parameters_gap-extend>
<Parameters_filter>L;m;</Parameters_filter>
</Parameters></BlastOutput_param>
<BlastOutput_iterations>
<Iteration>
<Iteration_iter-num>1</Iteration_iter-num>
<Iteration_query-ID>Query_1</Iteration_query-ID>
<Iteration_query-def>TestQuery</Iteration_query-def>
<Iteration_query-len>{query_len}</Iteration_query-len>
<Iteration_hits>
{hits_xml}
</Iteration_hits>
<Iteration_stat><Statistics>
<Statistics_db-num>70000000</Statistics_db-num>
<Statistics_db-len>100000000000</Statistics_db-len>
<Statistics_hsp-len>22</Statistics_hsp-len>
<Statistics_eff-space>9890000000000</Statistics_eff-space>
<Statistics_kappa>0.41</Statistics_kappa>
<Statistics_lambda>0.625</Statistics_lambda>
<Statistics_entropy>0.78</Statistics_entropy>
</Statistics></Iteration_stat>
</Iteration>
</BlastOutput_iterations>
</BlastOutput>""".encode()


def _no_hits_xml(query_len: int = 100) -> bytes:
    return f"""<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "NCBI_BlastOutput.dtd">
<BlastOutput>
<BlastOutput_program>blastn</BlastOutput_program>
<BlastOutput_version>BLASTN 2.14.0+</BlastOutput_version>
<BlastOutput_reference>Altschul SF et al.</BlastOutput_reference>
<BlastOutput_db>core_nt</BlastOutput_db>
<BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
<BlastOutput_query-def>TestQuery</BlastOutput_query-def>
<BlastOutput_query-len>{query_len}</BlastOutput_query-len>
<BlastOutput_param><Parameters>
<Parameters_expect>10</Parameters_expect>
<Parameters_sc-match>1</Parameters_sc-match>
<Parameters_sc-mismatch>-2</Parameters_sc-mismatch>
<Parameters_gap-open>0</Parameters_gap-open>
<Parameters_gap-extend>2</Parameters_gap-extend>
<Parameters_filter>L;m;</Parameters_filter>
</Parameters></BlastOutput_param>
<BlastOutput_iterations>
<Iteration>
<Iteration_iter-num>1</Iteration_iter-num>
<Iteration_query-ID>Query_1</Iteration_query-ID>
<Iteration_query-def>TestQuery</Iteration_query-def>
<Iteration_query-len>{query_len}</Iteration_query-len>
<Iteration_hits></Iteration_hits>
<Iteration_stat><Statistics>
<Statistics_db-num>70000000</Statistics_db-num>
<Statistics_db-len>100000000000</Statistics_db-len>
<Statistics_hsp-len>22</Statistics_hsp-len>
<Statistics_eff-space>9890000000000</Statistics_eff-space>
<Statistics_kappa>0.41</Statistics_kappa>
<Statistics_lambda>0.625</Statistics_lambda>
<Statistics_entropy>0.78</Statistics_entropy>
</Statistics></Iteration_stat>
</Iteration>
</BlastOutput_iterations>
</BlastOutput>""".encode()


def _parse_xml(xml_bytes: bytes):
    """Parse raw BLAST XML bytes into a Bio.Blast.Record."""
    from Bio import Blast
    return Blast.read(io.BytesIO(xml_bytes))


# ── Organism extraction ───────────────────────────────────────────────────────

class TestExtractOrganism:
    def test_standard_bracket_notation(self):
        assert _extract_organism("Some mRNA [Homo sapiens]") == "Homo sapiens"

    def test_trailing_whitespace(self):
        assert _extract_organism("Some mRNA [Mus musculus]  ") == "Mus musculus"

    def test_no_organism(self):
        assert _extract_organism("Some mRNA without organism") == ""

    def test_nested_brackets_picks_outermost_trailing(self):
        # Only the last [..] group at end of string should match
        result = _extract_organism("protein [isoform 1] mRNA [Arabidopsis thaliana]")
        assert result == "Arabidopsis thaliana"

    def test_empty_string(self):
        assert _extract_organism("") == ""


# ── parse_blast_record ────────────────────────────────────────────────────────

class TestParseBlastRecord:
    def test_basic_fields(self):
        record = _parse_xml(_blast_xml())
        hits = parse_blast_record(record)

        assert len(hits) == 1
        h = hits[0]
        assert h.accession == "NM_000001"
        assert "Homo sapiens" in h.description
        assert h.organism == "Homo sapiens"

    def test_query_length(self):
        record = _parse_xml(_blast_xml(query_len=100))
        hits = parse_blast_record(record)
        assert hits[0].query_length == 100

    def test_evalue(self):
        record = _parse_xml(_blast_xml(evalue="1.4e-45"))
        hits = parse_blast_record(record)
        assert hits[0].evalue == pytest.approx(1.4e-45, rel=1e-3)

    def test_bit_score(self):
        record = _parse_xml(_blast_xml(bit_score=180.5))
        hits = parse_blast_record(record)
        assert hits[0].bit_score == pytest.approx(180.5, rel=1e-3)

    def test_identity_count(self):
        record = _parse_xml(_blast_xml(identity=95, align_len=97))
        hits = parse_blast_record(record)
        assert hits[0].identity_count == 95

    def test_pct_identity(self):
        record = _parse_xml(_blast_xml(identity=95, align_len=97))
        hits = parse_blast_record(record)
        assert hits[0].pct_identity == pytest.approx(95 / 97, rel=1e-3)

    def test_query_coverage(self):
        # query is 100 bp; HSP covers positions 1-97 (0-based: 0-97)
        record = _parse_xml(_blast_xml(query_len=100, q_from=1, q_to=97, align_len=97))
        hits = parse_blast_record(record)
        assert hits[0].query_coverage == pytest.approx(0.97, rel=1e-2)

    def test_no_hits(self):
        record = _parse_xml(_no_hits_xml())
        hits = parse_blast_record(record)
        assert hits == []

    def test_multiple_hits(self):
        record = _parse_xml(_blast_xml(n_hits=3))
        hits = parse_blast_record(record)
        assert len(hits) == 3

    def test_alignment_text_not_empty(self):
        record = _parse_xml(_blast_xml())
        hits = parse_blast_record(record)
        assert len(hits[0].alignment_text) > 0

    def test_hsp_coordinates_stored(self):
        record = _parse_xml(_blast_xml(query_len=100, q_from=1, q_to=97, align_len=97))
        hits = parse_blast_record(record)
        h = hits[0]
        assert h.hsp_query_start == 0   # 0-based from Bio.Blast coords
        assert h.hsp_query_end == 97


# ── BLASTHit properties ───────────────────────────────────────────────────────

class TestBLASTHitProperties:
    def _make_hit(self, **kwargs) -> BLASTHit:
        defaults = dict(
            accession="NM_000001",
            description="desc [Homo sapiens]",
            organism="Homo sapiens",
            query_length=100,
            hsp_query_start=0,
            hsp_query_end=90,
            hsp_align_length=90,
            identity_count=85,
            evalue=1e-30,
            bit_score=150.0,
            score=80.0,
            alignment_text="...",
        )
        defaults.update(kwargs)
        return BLASTHit(**defaults)

    def test_query_coverage_normal(self):
        h = self._make_hit(query_length=100, hsp_query_start=0, hsp_query_end=90)
        assert h.query_coverage == pytest.approx(0.90)

    def test_query_coverage_zero_query_length(self):
        h = self._make_hit(query_length=0)
        assert h.query_coverage == 0.0

    def test_pct_identity_normal(self):
        h = self._make_hit(identity_count=85, hsp_align_length=90)
        assert h.pct_identity == pytest.approx(85 / 90)

    def test_pct_identity_zero_align_length(self):
        h = self._make_hit(hsp_align_length=0)
        assert h.pct_identity == 0.0

    def test_pct_identity_perfect(self):
        h = self._make_hit(identity_count=100, hsp_align_length=100)
        assert h.pct_identity == pytest.approx(1.0)


# ── MIN_QUERY_LENGTH constant ─────────────────────────────────────────────────

class TestConstants:
    def test_min_query_length_is_30(self):
        assert MIN_QUERY_LENGTH == 30

    def test_databases_has_core_nt_first(self):
        assert DATABASES[0][0] == "core_nt"

    def test_databases_has_three_entries(self):
        assert len(DATABASES) == 3

    def test_all_database_ids_are_strings(self):
        for db_id, _label in DATABASES:
            assert isinstance(db_id, str) and db_id


# ── Mocked BLASTWorker (no network) ──────────────────────────────────────────

class TestBLASTWorkerMocked:
    """Verify that BLASTWorker correctly uses parsed results and signals."""

    def _stringio_handle(self, xml_bytes: bytes) -> io.StringIO:
        """Return a StringIO — the exact type NCBIWWW.qblast() returns."""
        return io.StringIO(xml_bytes.decode("utf-8"))

    def test_worker_emits_results_on_success(self, qtbot=None):
        """Worker should parse XML and emit result_ready with hits list.

        The mock returns StringIO to match real NCBIWWW.qblast() behaviour
        (it does handle.read().decode() then wraps in StringIO).
        """
        try:
            from chromaview.gui.blast_dialog import BLASTWorker
        except ImportError:
            pytest.skip("PyQt6 not available")

        mock_handle = self._stringio_handle(_blast_xml())
        received = []

        worker = BLASTWorker(
            sequence="ATGCATGCATGCATGCATGCATGCATGCATGC",
            database="core_nt",
            email="test@example.com",
        )
        worker.result_ready.connect(received.append)
        errors = []
        worker.error.connect(errors.append)

        with patch("Bio.Blast.NCBIWWW.qblast", return_value=mock_handle):
            worker.run()

        assert not errors, f"Worker emitted error: {errors}"
        assert len(received) == 1
        hits = received[0]
        assert len(hits) == 1
        assert hits[0].accession == "NM_000001"

    def test_worker_emits_error_on_network_failure(self):
        try:
            from chromaview.gui.blast_dialog import BLASTWorker
        except ImportError:
            pytest.skip("PyQt6 not available")

        worker = BLASTWorker(
            sequence="ATGCATGC" * 10,
            database="core_nt",
            email="",
        )
        errors = []
        worker.error.connect(errors.append)

        with patch(
            "Bio.Blast.NCBIWWW.qblast",
            side_effect=OSError("Network unreachable"),
        ):
            worker.run()

        assert len(errors) == 1
        assert "Network unreachable" in errors[0]

    def test_worker_no_error_when_cancelled(self):
        try:
            from chromaview.gui.blast_dialog import BLASTWorker
        except ImportError:
            pytest.skip("PyQt6 not available")

        worker = BLASTWorker(
            sequence="ATGCATGC" * 10,
            database="core_nt",
            email="",
        )
        worker._cancelled = True
        errors = []
        worker.error.connect(errors.append)

        with patch(
            "Bio.Blast.NCBIWWW.qblast",
            side_effect=OSError("should be suppressed"),
        ):
            worker.run()

        assert not errors

    def test_worker_sets_ncbiwww_email_and_tool(self):
        try:
            from chromaview.gui.blast_dialog import BLASTWorker
            from Bio.Blast import NCBIWWW
        except ImportError:
            pytest.skip("PyQt6 not available")

        mock_handle = self._stringio_handle(_blast_xml())
        worker = BLASTWorker(
            sequence="ATGCATGC" * 10,
            database="nt",
            email="blast@example.com",
            tool="TestTool",
        )
        worker.result_ready.connect(lambda _: None)

        with patch("Bio.Blast.NCBIWWW.qblast", return_value=mock_handle):
            worker.run()

        assert NCBIWWW.email == "blast@example.com"
        assert NCBIWWW.tool == "TestTool"

    def test_empty_results_yields_empty_list(self):
        try:
            from chromaview.gui.blast_dialog import BLASTWorker
        except ImportError:
            pytest.skip("PyQt6 not available")

        mock_handle = self._stringio_handle(_no_hits_xml())
        received = []

        worker = BLASTWorker(
            sequence="ATGCATGC" * 10,
            database="core_nt",
            email="",
        )
        worker.result_ready.connect(received.append)

        with patch("Bio.Blast.NCBIWWW.qblast", return_value=mock_handle):
            worker.run()

        assert received == [[]]


# ── Mode-mismatch integration tests ──────────────────────────────────────────

class TestStreamModeHandling:
    """Verify the StringIO→BytesIO conversion that bridges qblast and Bio.Blast.read."""

    def test_stringio_raises_without_conversion(self):
        """Confirm StringIO is rejected by Bio.Blast.read — documents the root cause."""
        from Bio import Blast
        from Bio.Blast import StreamModeError

        xml_str = _blast_xml().decode("utf-8")
        with pytest.raises((StreamModeError, ValueError, Exception)):
            Blast.read(io.StringIO(xml_str))

    def test_stringio_works_after_encode_to_bytesio(self):
        """The worker's fix: encode StringIO content to BytesIO before parsing."""
        from Bio import Blast

        xml_str = _blast_xml().decode("utf-8")
        text_handle = io.StringIO(xml_str)          # what qblast actually returns
        raw = text_handle.read()
        if isinstance(raw, str):
            raw = raw.encode("utf-8")
        record = Blast.read(io.BytesIO(raw))        # must not raise
        hits = parse_blast_record(record)
        assert len(hits) == 1
        assert hits[0].accession == "NM_000001"

    def test_bytesio_still_works_directly(self):
        """BytesIO (e.g. from saved files) continues to parse without conversion."""
        from Bio import Blast

        record = Blast.read(io.BytesIO(_blast_xml()))
        hits = parse_blast_record(record)
        assert len(hits) == 1

    def test_conversion_preserves_all_hit_fields(self):
        """Full round-trip: StringIO → BytesIO conversion → parse → check fields.

        align_len, q_from, q_to, h_from, h_to must satisfy:
          q_to - q_from + 1 == align_len  (no gaps in our synthetic sequence)
        """
        from Bio import Blast

        xml_str = _blast_xml(
            query_len=150,
            hit_accession="AB123456",
            hit_def="Test organism gene [Arabidopsis thaliana]",
            evalue="3.2e-12",
            bit_score=95.3,
            identity=80,
            align_len=90,
            q_from=1,
            q_to=90,       # q_to - q_from + 1 == align_len
            h_from=101,
            h_to=190,      # h_to - h_from + 1 == align_len
            gaps=0,
        ).decode("utf-8")

        raw = io.StringIO(xml_str).read().encode("utf-8")
        record = Blast.read(io.BytesIO(raw))
        hits = parse_blast_record(record)

        assert hits[0].accession == "AB123456"
        assert hits[0].organism == "Arabidopsis thaliana"
        assert hits[0].evalue == pytest.approx(3.2e-12, rel=1e-3)
        assert hits[0].bit_score == pytest.approx(95.3, rel=1e-3)
        assert hits[0].identity_count == 80
        assert hits[0].query_length == 150
