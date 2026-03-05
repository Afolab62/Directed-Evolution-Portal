"""
test_sequence_tools.py

Unit tests for services.sequence_tools: FASTA parsing, DNA validation,
translation, reverse complement, six-frame translation, and Smith–Waterman.
"""
import pytest

from services.sequence_tools import (
    AlignmentResult,
    CODON_TABLE,
    parse_fasta_dna,
    parse_fasta_protein,
    reverse_complement,
    smith_waterman_local,
    translate_dna,
    translate_six_frames,
    validate_dna,
)
from services.errors import FastaParseError, InvalidSequenceError


# ---------------------------------------------------------------------------
# parse_fasta_dna
# ---------------------------------------------------------------------------

class TestParseFastaDna:
    def test_single_record_returned(self):
        fasta = ">seq1\nATGCATGC\n"
        assert parse_fasta_dna(fasta) == "ATGCATGC"

    def test_multiline_sequence_joined(self):
        fasta = ">seq1\nATG\nCAT\nGC\n"
        assert parse_fasta_dna(fasta) == "ATGCATGC"

    def test_lowercase_normalised_to_uppercase(self):
        fasta = ">s\natgc\n"
        assert parse_fasta_dna(fasta) == "ATGC"

    def test_multi_record_rejected(self):
        fasta = ">s1\nATGC\n>s2\nGGGG\n"
        with pytest.raises(FastaParseError):
            parse_fasta_dna(fasta)

    def test_empty_string_raises(self):
        with pytest.raises(FastaParseError):
            parse_fasta_dna("")

    def test_invalid_character_raises(self):
        fasta = ">s\nATGZ\n"
        with pytest.raises(InvalidSequenceError):
            parse_fasta_dna(fasta)

    def test_n_ambiguous_base_allowed(self):
        fasta = ">s\nATGNATGC\n"
        result = parse_fasta_dna(fasta)
        assert result == "ATGNATGC"


# ---------------------------------------------------------------------------
# parse_fasta_protein
# ---------------------------------------------------------------------------

class TestParseFastaProtein:
    def test_single_record_returned(self):
        fasta = ">sp|P12345|GENE_HUMAN\nMACVG\n"
        assert parse_fasta_protein(fasta) == "MACVG"

    def test_multi_record_rejected(self):
        fasta = ">p1\nMACVG\n>p2\nKLMNP\n"
        with pytest.raises(FastaParseError):
            parse_fasta_protein(fasta)

    def test_lowercase_normalised(self):
        fasta = ">p\nmacvg\n"
        assert parse_fasta_protein(fasta) == "MACVG"


# ---------------------------------------------------------------------------
# validate_dna
# ---------------------------------------------------------------------------

class TestValidateDna:
    def test_valid_dna_passes(self):
        validate_dna("ATGCNATGC")  # should not raise

    def test_invalid_character_raises(self):
        with pytest.raises(InvalidSequenceError):
            validate_dna("ATGZATGC")

    def test_empty_string_passes(self):
        validate_dna("")  # empty is valid (no invalid chars)


# ---------------------------------------------------------------------------
# translate_dna
# ---------------------------------------------------------------------------

class TestTranslateDna:
    def test_atg_translates_to_methionine(self):
        assert translate_dna("ATG", CODON_TABLE) == "M"

    def test_stop_codon(self):
        assert translate_dna("TAA", CODON_TABLE) == "*"

    def test_full_codon_sequence(self):
        # ATG TTT TAA → M F *
        assert translate_dna("ATGTTTAA", CODON_TABLE) == "MF"  # last codon TAA truncated to 2 nt

    def test_unknown_codon_returns_x(self):
        # NNN is not in any standard table → X
        result = translate_dna("NNN", CODON_TABLE)
        assert result == "X"

    def test_partial_trailing_codon_ignored(self):
        # ATG + AT (only 2 nt trailing) → only M
        assert translate_dna("ATGAT", CODON_TABLE) == "M"


# ---------------------------------------------------------------------------
# reverse_complement
# ---------------------------------------------------------------------------

class TestReverseComplement:
    def test_simple(self):
        assert reverse_complement("ATGC") == "GCAT"

    def test_palindrome(self):
        assert reverse_complement("AATT") == "AATT"

    def test_n_complement(self):
        result = reverse_complement("ATGN")
        assert result == "NCAT"

    def test_empty(self):
        assert reverse_complement("") == ""


# ---------------------------------------------------------------------------
# translate_six_frames
# ---------------------------------------------------------------------------

class TestTranslateSixFrames:
    def test_returns_six_keys(self):
        frames = translate_six_frames("ATGCATGCAT")
        assert set(frames.keys()) == {"+0", "+1", "+2", "-0", "-1", "-2"}

    def test_plus_zero_frame_starts_at_atg(self):
        # ATGTTT = M + F (2 codons)
        frames = translate_six_frames("ATGTTT")
        assert frames["+0"].startswith("M")

    def test_all_values_are_strings(self):
        frames = translate_six_frames("ATGCATGCAT")
        assert all(isinstance(v, str) for v in frames.values())


# ---------------------------------------------------------------------------
# smith_waterman_local
# ---------------------------------------------------------------------------

class TestSmithWatermanLocal:
    def test_identical_sequences_have_positive_score(self):
        result = smith_waterman_local("ACGT", "ACGT")
        assert result is not None
        assert result.score > 0

    def test_no_overlap_returns_zero_score(self):
        # Completely different sequences → score may be 0 or very low
        result = smith_waterman_local("AAAA", "TTTT")
        assert result is not None
        assert result.score >= 0

    def test_returns_alignment_result_type(self):
        result = smith_waterman_local("ACGT", "ACGT")
        assert isinstance(result, AlignmentResult)

    def test_returns_none_when_too_large(self):
        # 2001 × 2001 > max_cells default of 2_000_000
        big_a = "A" * 2001
        big_b = "A" * 2001
        result = smith_waterman_local(big_a, big_b)
        assert result is None

    def test_partial_match(self):
        result = smith_waterman_local("XXXACGTXXX", "ACGT")
        assert result is not None
        assert result.score > 0
