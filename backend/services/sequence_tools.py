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


# =============================================================================
# FASTA parsing
# =============================================================================

def _parse_fasta(text: str) -> List[str]:
    # strips blank lines and header lines starting with '>' before building
    # sequences — handles FASTA files with inconsistent line endings
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if not lines:
        raise FastaParseError("Empty FASTA input")

    sequences: List[str] = []
    current: List[str] = []

    for line in lines:
        if line.startswith(">"):
            # on a new header, flush the accumulated sequence before starting the next
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
    # enforces single-record input — multi-record FASTA would be ambiguous
    # for plasmid validation, which expects exactly one reference sequence
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
    # 'N' is included as a valid character — it represents an ambiguous base
    # commonly produced by sequencing when the instrument cannot confidently
    # call a nucleotide; rejecting N would fail many real sequencing outputs
    invalid = set(seq.upper()) - {"A", "C", "G", "T", "N"}
    if invalid:
        raise InvalidSequenceError(
            f"Invalid DNA characters found: {', '.join(sorted(invalid))}"
        )


# =============================================================================
# Translation
# =============================================================================

# standard genetic code mapping all 64 codons to single-letter amino acids;
# TAA, TAG, TGA are stop codons represented as '*'
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
    # range stops at len - 2 to avoid reading a partial codon at the end of
    # the sequence — a trailing 1 or 2 nucleotide remainder is silently ignored
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        # unknown codons (e.g. those containing 'N') are translated as 'X'
        # rather than raising an error, consistent with standard bioinformatics convention
        aa.append(codon_table.get(codon, "X"))
    return "".join(aa)


def reverse_complement(seq: str) -> str:
    # reads the sequence in reverse then maps each base to its Watson-Crick
    # complement — required to search the antisense strand of a DNA molecule
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(complement.get(b, "N") for b in reversed(seq))


def translate_six_frames(
    seq: str,
    *,
    codon_table: Dict[str, str] | None = None,
) -> Dict[str, str]:
    # a double-stranded DNA molecule has 6 possible reading frames:
    # 3 on the sense strand (+0, +1, +2) and 3 on the antisense strand (-0, -1, -2).
    # genes can be encoded on either strand, so all 6 must be searched
    if codon_table is None:
        codon_table = CODON_TABLE

    frames: Dict[str, str] = {}

    for frame in range(3):
        frames[f"+{frame}"] = translate_dna(seq[frame:], codon_table)

    rev = reverse_complement(seq)
    for frame in range(3):
        frames[f"-{frame}"] = translate_dna(rev[frame:], codon_table)

    return frames


# =============================================================================
# Smith–Waterman local alignment
# =============================================================================

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
    # local alignment is O(n*m) in both time and memory — for very long
    # sequences (e.g. whole plasmids vs large proteins) this becomes prohibitive.
    # returning None signals to the caller that alignment was skipped
    if n * m > max_cells:
        return None

    # H[i][j] stores the best local alignment score ending at position (i, j);
    # initialised to zero so the alignment can start anywhere in either sequence
    H = [[0] * (m + 1) for _ in range(n + 1)]

    best = AlignmentResult(0, 0, 0)

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score = max(
                0,  # zero floor allows the alignment to restart — key property of local alignment
                H[i - 1][j - 1] + (match if a[i - 1] == b[j - 1] else mismatch),
                H[i - 1][j] + gap,   # gap in b
                H[i][j - 1] + gap,   # gap in a
            )
            H[i][j] = score
            if score > best.score:
                best = AlignmentResult(score, i, j)

    return best
