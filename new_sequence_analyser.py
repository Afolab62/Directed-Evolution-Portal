import sqlite3
from typing import List, Tuple, Optional
from dataclasses import dataclass
import logging

from Bio import Align  #Biopython PairwiseAligner

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


@dataclass
class MutationInfo:
    """Data class to store mutation information"""
    position: int
    wt_codon: str
    mut_codon: str
    wt_aa: str
    mut_aa: str
    mutation_type: str  # 'synonymous' or 'non-synonymous'

    def __str__(self):
        return f"{self.wt_aa}{self.position}{self.mut_aa} ({self.mutation_type})"


class SequenceAnalyzer:
    """Sequence analysis and mutation detection using Biopython Smith-Waterman local alignment"""

    # Standard genetic code
    GENETIC_CODE = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

    def __init__(self, wt_protein_seq: str, wt_plasmid_seq: str):
        """Initialize with wild-type sequences"""
        self.wt_protein_seq = wt_protein_seq.upper().strip()
        self.wt_plasmid_seq = wt_plasmid_seq.upper().strip()

        self.wt_gene_seq: Optional[str] = None
        self.wt_gene_start: Optional[int] = None  # 0-based start, on the provided plasmid sequence
        self.wt_gene_strand: Optional[str] = None  # '+' or '-'

        # Configure a reusable Biopython aligner (local alignment).
        # With a linear gap penalty, PairwiseAligner will select Smith-Waterman for local alignment.
        self._aligner = Align.PairwiseAligner()
        self._aligner.mode = "local"
        self._aligner.match_score = 2.0
        self._aligner.mismatch_score = -1.0
        self._aligner.gap_score = -1.0

        self._locate_wt_gene()

    def _smith_waterman(self, seq1: str, seq2: str, match=2, mismatch=-1, gap=-1) -> Tuple[float, int, int]:
        """
        Smith-Waterman local alignment using Biopython PairwiseAligner.

        Returns: (score, start_pos_in_seq1, end_pos_in_seq1)
        where start/end are 0-based indices in seq1, end-exclusive.
        """
        # Keep configuration explicit to avoid version-dependent defaults.
        if (self._aligner.mode != "local"
                or self._aligner.match_score != float(match)
                or self._aligner.mismatch_score != float(mismatch)
                or self._aligner.gap_score != float(gap)):
            self._aligner.mode = "local"
            self._aligner.match_score = float(match)
            self._aligner.mismatch_score = float(mismatch)
            self._aligner.gap_score = float(gap)

        # Optional sanity check: confirm Biopython selected Smith-Waterman.
        # (Do not hard-fail: algorithm choice can depend on gap scoring configuration.)
        if getattr(self._aligner, "algorithm", None) != "Smith-Waterman":
            logger.debug(f"PairwiseAligner selected algorithm='{self._aligner.algorithm}'")

        alignments = self._aligner.align(seq1, seq2)
        if len(alignments) == 0:
            return 0.0, 0, 0

        alignment = alignments[0]  # best-scoring alignment (one of them, if multiple optima)

        # alignment.aligned gives the aligned blocks for target (row 0) and query (row 1).
        target_blocks = alignment.aligned[0]
        start_pos = int(target_blocks[0][0])
        end_pos = int(target_blocks[-1][1])

        return float(alignment.score), start_pos, end_pos

    def _translate(self, dna_seq: str) -> str:
        """Translate DNA to protein sequence"""
        protein = []
        for i in range(0, len(dna_seq) - 2, 3):
            codon = dna_seq[i:i+3]
            aa = self.GENETIC_CODE.get(codon, 'X')
            if aa == '*':
                break
            protein.append(aa)
        return ''.join(protein)

    def _reverse_complement(self, seq: str) -> str:
        """Generate reverse complement of DNA sequence"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return ''.join(complement.get(base, base) for base in reversed(seq))

    def _locate_wt_gene(self):
        """Locate gene by local alignment on translated sequences (both strands, all frames)"""
        logger.info("Locating WT gene using Biopython PairwiseAligner (Smith-Waterman local alignment)...")

        plasmid_len = len(self.wt_plasmid_seq)
        extended_plasmid = self.wt_plasmid_seq + self.wt_plasmid_seq  # circularisation trick
        ext_len = len(extended_plasmid)

        best_score = float("-inf")
        best_match = None

        # Try both strands
        for strand, strand_seq in [("+", extended_plasmid), ("-", self._reverse_complement(extended_plasmid))]:
            # Try all three reading frames
            for frame in range(3):
                translated = self._translate(strand_seq[frame:])

                score, aa_start, aa_end = self._smith_waterman(translated, self.wt_protein_seq)

                if score > best_score:
                    dna_start = frame + (aa_start * 3)
                    dna_end = frame + (aa_end * 3)
                    gene_seq = strand_seq[dna_start:dna_end]

                    if strand == "+":
                        gene_start_on_input = dna_start % plasmid_len
                    else:
                        # Map reverse-complement coordinates back to the original (input) plasmid.
                        # If rc[i] = complement(ext[ext_len-1-i]), then:
                        # rc[a:b] corresponds to ext[ext_len-b : ext_len-a] on the input strand.
                        ext_start = ext_len - dna_end
                        gene_start_on_input = ext_start % plasmid_len

                    best_score = score
                    best_match = {
                        "score": score,
                        "strand": strand,
                        "frame": frame,
                        "gene_seq": gene_seq,
                        "start_on_input": gene_start_on_input,
                        "protein": self._translate(gene_seq)
                    }

        # Threshold borrowed from the original logic; with match=2 a perfect match is 2*L
        if not best_match or best_match["score"] < len(self.wt_protein_seq) * 1.5:
            raise ValueError("Could not locate WT gene in plasmid!")

        self.wt_gene_seq = best_match["gene_seq"]
        self.wt_gene_start = best_match["start_on_input"]
        self.wt_gene_strand = best_match["strand"]

        logger.info(
            f"WT gene found: score={best_match['score']}, start={self.wt_gene_start}, "
            f"strand={self.wt_gene_strand}, length={len(self.wt_gene_seq)}"
        )

    def extract_gene_from_variant(self, variant_plasmid_seq: str) -> str:
        """Extract gene from variant plasmid using coordinates and strand determined from WT"""
        variant_plasmid_seq = variant_plasmid_seq.upper().strip()
        gene_len = len(self.wt_gene_seq)

        # Extract the segment from the INPUT strand coordinates (handle circular DNA)
        if self.wt_gene_start + gene_len <= len(variant_plasmid_seq):
            segment = variant_plasmid_seq[self.wt_gene_start:self.wt_gene_start + gene_len]
        else:
            segment = (
                variant_plasmid_seq[self.wt_gene_start:]
                + variant_plasmid_seq[:gene_len - (len(variant_plasmid_seq) - self.wt_gene_start)]
            )

        # If the gene is on the reverse strand, return it in coding orientation
        if self.wt_gene_strand == "-":
            segment = self._reverse_complement(segment)

        return segment

    def identify_mutations(self, variant_gene: str) -> List[MutationInfo]:
        """Identify mutations comparing variant to WT"""
        mutations: List[MutationInfo] = []

        if len(variant_gene) != len(self.wt_gene_seq):
            logger.warning(f"Length mismatch: variant={len(variant_gene)}, WT={len(self.wt_gene_seq)}")
            return mutations

        for i in range(0, len(self.wt_gene_seq) - 2, 3):
            wt_codon = self.wt_gene_seq[i:i+3]
            mut_codon = variant_gene[i:i+3]

            if wt_codon != mut_codon:
                wt_aa = self.GENETIC_CODE.get(wt_codon, 'X')
                mut_aa = self.GENETIC_CODE.get(mut_codon, 'X')

                mutations.append(
                    MutationInfo(
                        position=(i // 3) + 1,
                        wt_codon=wt_codon,
                        mut_codon=mut_codon,
                        wt_aa=wt_aa,
                        mut_aa=mut_aa,
                        mutation_type="synonymous" if wt_aa == mut_aa else "non-synonymous",
                    )
                )

        return mutations

    def analyze_variant(self, variant_plasmid: str) -> dict:
        """Complete analysis of a variant plasmid"""
        gene_seq = self.extract_gene_from_variant(variant_plasmid)
        protein_seq = self._translate(gene_seq)
        mutations = self.identify_mutations(gene_seq)

        num_synonymous = sum(1 for m in mutations if m.mutation_type == 'synonymous')
        num_nonsynonymous = len(mutations) - num_synonymous

        logger.info(
            f"Analysis complete: {len(mutations)} mutations "
            f"({num_synonymous} syn, {num_nonsynonymous} non-syn)"
        )

        return {
            'gene_sequence': gene_seq,
            'protein_sequence': protein_seq,
            'mutations': mutations,
            'num_mutations': len(mutations),
            'num_synonymous': num_synonymous,
            'num_nonsynonymous': num_nonsynonymous
        }
