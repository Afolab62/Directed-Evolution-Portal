"""
Sequence utility functions used by plasmid validation.

Includes:
- FASTA parsing
- DNA validation
- Translation utilities
- Smith–Waterman local alignment
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List

from app.services.errors import FastaParseError, InvalidSequenceError


# =========================
# FASTA parsing
# =========================

def _parse_fasta(text: str) -> List[str]:
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if not lines:
        raise FastaParseError("Empty FASTA input")

    sequences: List[str] = []
    current: List[str] = []

    for line in lines:
        if line.startswith(">"):
            if current:
                sequences.append("".join(current))
                current = []
        else:
            current.append(line)

    if current:
        sequences.append("".join(current))

    if not sequences:
        raise FastaParseError("No sequences found in FASTA")

    return sequences


def parse_fasta_dna(text: str) -> str:
    sequences = _parse_fasta(text)
    if len(sequences) != 1:
        raise FastaParseError("Expected exactly one DNA sequence")

    seq = sequences[0].upper()
    validate_dna(seq)
    return seq


def parse_fasta_protein(text: str) -> str:
    sequences = _parse_fasta(text)
    if len(sequences) != 1:
        raise FastaParseError("Expected exactly one protein sequence")

    return sequences[0].upper()


def validate_dna(seq: str) -> None:
    invalid = set(seq.upper()) - {"A", "C", "G", "T", "N"}
    if invalid:
        raise InvalidSequenceError(
            f"Invalid DNA characters found: {', '.join(sorted(invalid))}"
        )


# =========================
# Translation
# =========================

CODON_TABLE: Dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate_dna(seq: str, codon_table: Dict[str, str]) -> str:
    aa: List[str] = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        aa.append(codon_table.get(codon, "X"))
    return "".join(aa)


def reverse_complement(seq: str) -> str:
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(complement.get(b, "N") for b in reversed(seq))


def translate_six_frames(
    seq: str,
    *,
    codon_table: Dict[str, str] | None = None,
) -> Dict[str, str]:
    if codon_table is None:
        codon_table = CODON_TABLE

    frames: Dict[str, str] = {}

    for frame in range(3):
        frames[f"+{frame}"] = translate_dna(seq[frame:], codon_table)

    rev = reverse_complement(seq)
    for frame in range(3):
        frames[f"-{frame}"] = translate_dna(rev[frame:], codon_table)

    return frames


# =========================
# Smith–Waterman alignment
# =========================

@dataclass
class AlignmentResult:
    score: int
    end_i: int
    end_j: int


def smith_waterman_local(
    a: str,
    b: str,
    *,
    match: int = 2,
    mismatch: int = -1,
    gap: int = -1,
    max_cells: int = 2_000_000,
) -> AlignmentResult | None:
    n, m = len(a), len(b)
    if n * m > max_cells:
        return None

    H = [[0] * (m + 1) for _ in range(n + 1)]

    best = AlignmentResult(0, 0, 0)

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score = max(
                0,
                H[i - 1][j - 1] + (match if a[i - 1] == b[j - 1] else mismatch),
                H[i - 1][j] + gap,
                H[i][j - 1] + gap,
            )
            H[i][j] = score
            if score > best.score:
                best = AlignmentResult(score, i, j)

    return best
