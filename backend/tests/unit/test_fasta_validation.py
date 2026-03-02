import pytest

from app.services.sequence_tools import parse_fasta_dna, parse_fasta_protein
from app.services.errors import FastaParseError, InvalidSequenceError


def test_multi_record_dna_fasta_rejected_by_default():
    fasta = """>seq1
ACGTACGT
>seq2
ACGTACGT
"""
    with pytest.raises(FastaParseError):
        parse_fasta_dna(fasta)


def test_invalid_dna_character_rejected():
    fasta = """>plasmid
ACGTZACGT
"""
    with pytest.raises(InvalidSequenceError):
        parse_fasta_dna(fasta)


def test_multi_record_protein_fasta_rejected_by_default():
    fasta = """>p1
ACDEFG
>p2
ACDEFG
"""
    with pytest.raises(FastaParseError):
        parse_fasta_protein(fasta)
