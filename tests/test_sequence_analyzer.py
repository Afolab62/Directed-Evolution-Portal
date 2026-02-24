"""
tests/test_sequence_analyzer.py
---------------------------------
Unit tests for SequenceAnalyzer.

Key point: because SequenceAnalyzer accepts a StagingResult object
rather than fetching it from the DB itself, ALL tests here run with
zero database involvement — just plain Python objects.

Run with:  pytest tests/test_sequence_analyzer.py -v
"""

import pytest
from database.models import StagingResult, MutationInfo
from analysis.sequence_analyzer import SequenceAnalyzer


# ---------------------------------------------------------------------------
# Minimal fixtures
# ---------------------------------------------------------------------------

# Short synthetic WT gene: encodes M-A-K-* (Met-Ala-Lys-stop)
WT_GENE = "ATGGCTAAATAA"          # 12 bp, 3 coding codons
WT_PROTEIN = "MAK"

# Pretend the gene sits at position 100 in a 200 bp plasmid
PLASMID_PREFIX = "N" * 100
PLASMID_SUFFIX = "N" * (200 - 100 - len(WT_GENE))
WT_PLASMID = PLASMID_PREFIX + WT_GENE + PLASMID_SUFFIX


def make_staging_result(
    start: int = 100,
    end: int = 112,
    strand: str = "+",
    frame: int = 0,
    wraps_origin: bool = False,
) -> StagingResult:
    """Helper: build a minimal valid StagingResult for testing."""
    return StagingResult(
        id=1,
        experiment_id=1,
        is_valid=True,
        match_type="exact",
        identity=1.0,
        coverage=1.0,
        strand=strand,
        frame=frame,
        start_nt=start,
        end_nt_exclusive=end,
        wraps_origin=wraps_origin,
        notes=None,
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestSequenceAnalyzerWithStagingResult:
    """All tests inject a StagingResult — no DB required."""

    def test_loads_wt_gene_correctly(self):
        sr = make_staging_result()
        analyzer = SequenceAnalyzer(WT_PROTEIN, WT_PLASMID, staging_result=sr)
        assert analyzer.wt_gene_seq == WT_GENE

    def test_no_mutations_on_wt(self):
        sr = make_staging_result()
        analyzer = SequenceAnalyzer(WT_PROTEIN, WT_PLASMID, staging_result=sr)
        result = analyzer.analyze_variant(WT_PLASMID)
        assert result["num_mutations"] == 0
        assert result["mutations"] == []

    def test_detects_nonsynonymous_mutation(self):
        # Change codon 2: GCT (Ala) → GTT (Val) — non-synonymous
        mutant_gene = "ATGGTTAAATAA"
        mutant_plasmid = PLASMID_PREFIX + mutant_gene + PLASMID_SUFFIX

        sr = make_staging_result()
        analyzer = SequenceAnalyzer(WT_PROTEIN, WT_PLASMID, staging_result=sr)
        result = analyzer.analyze_variant(mutant_plasmid)

        assert result["num_mutations"] == 1
        assert result["num_nonsynonymous"] == 1
        assert result["num_synonymous"] == 0

        mut: MutationInfo = result["mutations"][0]
        assert mut.position == 2
        assert mut.wt_aa == "A"
        assert mut.mut_aa == "V"
        assert mut.mutation_type == "non-synonymous"

    def test_detects_synonymous_mutation(self):
        # Change codon 2: GCT → GCC — both encode Ala (synonymous)
        mutant_gene = "ATGGCCAAATAA"
        mutant_plasmid = PLASMID_PREFIX + mutant_gene + PLASMID_SUFFIX

        sr = make_staging_result()
        analyzer = SequenceAnalyzer(WT_PROTEIN, WT_PLASMID, staging_result=sr)
        result = analyzer.analyze_variant(mutant_plasmid)

        assert result["num_mutations"] == 1
        assert result["num_synonymous"] == 1
        assert result["num_nonsynonymous"] == 0
        assert result["mutations"][0].mutation_type == "synonymous"

    def test_raises_on_invalid_staging_result(self):
        sr = make_staging_result()
        sr.is_valid = False
        sr.notes = "alignment failed"
        with pytest.raises(ValueError, match="is_valid=False"):
            SequenceAnalyzer(WT_PROTEIN, WT_PLASMID, staging_result=sr)

    def test_raises_on_length_mismatch(self):
        """Indel in variant gene should raise, not silently produce garbage."""
        sr = make_staging_result()
        analyzer = SequenceAnalyzer(WT_PROTEIN, WT_PLASMID, staging_result=sr)

        # Insert one extra base → length mismatch
        indel_gene = "ATGGCTAAATAAC"
        indel_plasmid = PLASMID_PREFIX + indel_gene + PLASMID_SUFFIX[:len(PLASMID_SUFFIX)-1]

        with pytest.raises(ValueError, match="length mismatch"):
            analyzer.identify_mutations(indel_gene)
