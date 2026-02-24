"""
analysis/sequence_analyzer.py
------------------------------
Pure sequence analysis and mutation detection.

Design principle
----------------
This class does ONE thing: take sequences in, return analysis results.
It has NO knowledge of the database, Flask, or any other infrastructure.
The staging result is injected by the service layer (app/services/analysis_service.py),
keeping this class fully unit-testable without a database.

If no StagingResult is injected, the class falls back to a Smith-Waterman
search — slower but always available as a safety net.
"""

import logging
from typing import List, Tuple, Optional

from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

from database.models import MutationInfo, StagingResult


logger = logging.getLogger(__name__)


class SequenceAnalyzer:
    """
    Sequence analysis and mutation detection.

    Gene coordinates are injected via a StagingResult from the plasmid
    validation pipeline.  _locate_wt_gene() is only called as a fallback
    when no staging result is provided.
    """

    # Standard genetic code (codon → single-letter AA)
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
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    }

    def __init__(
        self,
        wt_protein_seq: str,
        wt_plasmid_seq: str,
        staging_result: Optional[StagingResult] = None,
    ):
        """
        Parameters
        ----------
        wt_protein_seq : str
            Wild-type protein sequence (single-letter AA codes).
        wt_plasmid_seq : str
            Wild-type plasmid sequence (full circular DNA).
        staging_result : StagingResult, optional
            Pre-validated gene location from the plasmid validation pipeline.
            Provided by the service layer — this class never fetches it itself.
            When None, falls back to Smith-Waterman alignment (slower).
        """
        self.wt_protein_seq = wt_protein_seq.upper().strip()
        self.wt_plasmid_seq = wt_plasmid_seq.upper().strip()
        self.wt_gene_seq: Optional[str] = None
        self.wt_gene_start: Optional[int] = None
        self.staging_result = staging_result

        if staging_result is not None:
            self._load_from_staging(staging_result)
        else:
            logger.warning(
                "No StagingResult provided — falling back to Smith-Waterman. "
                "Consider running the plasmid validation pipeline first for speed and reliability."
            )
            self._locate_wt_gene()

    # ------------------------------------------------------------------
    # Staging result integration
    # ------------------------------------------------------------------

    def _load_from_staging(self, sr: StagingResult) -> None:
        """
        Populate wt_gene_seq / wt_gene_start directly from a StagingResult.
        Raises ValueError if the staging result is marked invalid.
        """
        if not sr.is_valid:
            raise ValueError(
                f"StagingResult id={sr.id} is marked is_valid=False "
                f"(notes: {sr.notes}). Cannot proceed — re-run plasmid validation."
            )

        plasmid = self.wt_plasmid_seq
        start = sr.start_nt
        end = sr.end_nt_exclusive

        # Handle genes that wrap around the plasmid origin
        if sr.wraps_origin:
            gene_seq = plasmid[start:] + plasmid[:end % len(plasmid)]
        else:
            gene_seq = plasmid[start:end]

        # Reverse-complement for minus-strand genes
        if sr.strand == '-':
            gene_seq = self._reverse_complement(gene_seq)

        # Correct reading frame then trim to codon boundary
        gene_seq = gene_seq[sr.frame:]
        gene_seq = gene_seq[:len(gene_seq) - len(gene_seq) % 3]

        self.wt_gene_seq = gene_seq
        self.wt_gene_start = start

        logger.info(
            f"Gene loaded from StagingResult id={sr.id}: "
            f"strand={sr.strand}, frame={sr.frame}, start={start}, "
            f"length={len(gene_seq)} bp, identity={sr.identity:.3f}, "
            f"coverage={sr.coverage:.3f}"
        )

    # ------------------------------------------------------------------
    # Fallback alignment helpers (used only when no StagingResult given)
    # ------------------------------------------------------------------

    def _smith_waterman(
        self,
        seq1: str,
        seq2: str,
        match: int = 2,
        mismatch: int = -1,
        gap: int = -1,
    ) -> Tuple[int, int, int]:
        """
        Smith-Waterman local alignment via BioPython PairwiseAligner.
        Returns (score, start_in_seq1, end_in_seq1).
        """
        aligner = PairwiseAligner()
        aligner.mode = "local"
        aligner.match_score = match
        aligner.mismatch_score = mismatch
        aligner.open_gap_score = gap
        aligner.extend_gap_score = gap

        alignments = aligner.align(seq1, seq2)
        if len(alignments) == 0:
            return 0, 0, 0

        best = alignments[0]
        target_blocks = best.aligned[0]
        if len(target_blocks) == 0:
            return int(best.score), 0, 0

        start = int(target_blocks[0][0])
        end = int(target_blocks[-1][1])
        return int(best.score), start, end

    def _translate(self, dna_seq: str) -> str:
        """Translate a DNA sequence to protein using BioPython (stops at first stop codon)."""
        trimmed = dna_seq[:len(dna_seq) - len(dna_seq) % 3]
        return str(Seq(trimmed).translate(to_stop=True))

    def _reverse_complement(self, seq: str) -> str:
        """Return the reverse complement of a DNA sequence using BioPython."""
        return str(Seq(seq).reverse_complement())

    def _locate_wt_gene(self) -> None:
        """
        Fallback: locate the gene via Smith-Waterman on all 6 reading frames.
        Only called when no StagingResult is provided.
        Handles circular plasmids by doubling the sequence before searching.
        """
        logger.info("Locating WT gene via Smith-Waterman (6-frame search)...")

        # Double the sequence to handle genes that wrap the origin
        extended_plasmid = self.wt_plasmid_seq + self.wt_plasmid_seq
        best_score = 0
        best_match = None

        for strand_seq in [extended_plasmid, self._reverse_complement(extended_plasmid)]:
            for frame in range(3):
                translated = self._translate(strand_seq[frame:])
                score, start, end = self._smith_waterman(translated, self.wt_protein_seq)

                if score > best_score:
                    best_score = score
                    dna_start = frame + (start * 3)
                    dna_end = frame + (end * 3)
                    gene_seq = strand_seq[dna_start:dna_end]
                    best_match = {
                        'score': score,
                        'gene_seq': gene_seq,
                        'start': dna_start % len(self.wt_plasmid_seq),
                    }

        if not best_match or best_match['score'] < len(self.wt_protein_seq) * 1.5:
            raise ValueError(
                "Could not locate WT gene in plasmid via Smith-Waterman. "
                "Check that the correct plasmid and protein were provided."
            )

        self.wt_gene_seq = best_match['gene_seq']
        self.wt_gene_start = best_match['start']

        logger.info(
            f"WT gene located: score={best_match['score']}, "
            f"start={self.wt_gene_start}, length={len(self.wt_gene_seq)} bp"
        )

    # ------------------------------------------------------------------
    # Core analysis — these are the public methods
    # ------------------------------------------------------------------

    def extract_gene_from_variant(self, variant_plasmid_seq: str) -> str:
        """
        Extract the gene region from a variant plasmid using WT coordinates.
        Handles origin-wrapping and strand/frame correction from staging.
        """
        variant_plasmid_seq = variant_plasmid_seq.upper().strip()
        gene_len = len(self.wt_gene_seq)

        # Determine if the gene wraps the plasmid origin
        wraps = (
            self.staging_result.wraps_origin
            if self.staging_result is not None
            else (self.wt_gene_start + gene_len > len(variant_plasmid_seq))
        )

        if not wraps:
            gene_seq = variant_plasmid_seq[self.wt_gene_start:self.wt_gene_start + gene_len]
        else:
            overflow = gene_len - (len(variant_plasmid_seq) - self.wt_gene_start)
            gene_seq = variant_plasmid_seq[self.wt_gene_start:] + variant_plasmid_seq[:overflow]

        # Apply strand / frame correction from staging result if available
        if self.staging_result is not None:
            if self.staging_result.strand == '-':
                gene_seq = self._reverse_complement(gene_seq)
            gene_seq = gene_seq[self.staging_result.frame:]
            gene_seq = gene_seq[:len(gene_seq) - len(gene_seq) % 3]

        return gene_seq

    def identify_mutations(self, variant_gene: str) -> List[MutationInfo]:
        """
        Identify mutations by comparing variant gene to WT codon by codon.

        Raises ValueError if the sequences differ in length
        (indicates an indel or frameshift — codon comparison cannot proceed).
        """
        mutations: List[MutationInfo] = []

        if len(variant_gene) != len(self.wt_gene_seq):
            raise ValueError(
                f"Gene length mismatch (variant={len(variant_gene)} bp, "
                f"WT={len(self.wt_gene_seq)} bp). "
                "This likely indicates an indel or frameshift. "
                "Codon-level comparison cannot proceed."
            )

        for i in range(0, len(self.wt_gene_seq) - 2, 3):
            wt_codon = self.wt_gene_seq[i:i + 3]
            mut_codon = variant_gene[i:i + 3]

            if wt_codon != mut_codon:
                wt_aa = self.GENETIC_CODE.get(wt_codon, 'X')
                mut_aa = self.GENETIC_CODE.get(mut_codon, 'X')
                mutations.append(MutationInfo(
                    position=(i // 3) + 1,   # 1-based
                    wt_codon=wt_codon,
                    mut_codon=mut_codon,
                    wt_aa=wt_aa,
                    mut_aa=mut_aa,
                    mutation_type='synonymous' if wt_aa == mut_aa else 'non-synonymous',
                ))

        return mutations

    def analyze_variant(self, variant_plasmid: str) -> dict:
        """
        Full analysis of a single variant plasmid sequence.

        Returns
        -------
        dict with keys:
            gene_sequence, protein_sequence, mutations (List[MutationInfo]),
            num_mutations, num_synonymous, num_nonsynonymous
        """
        gene_seq = self.extract_gene_from_variant(variant_plasmid)
        protein_seq = self._translate(gene_seq)
        mutations = self.identify_mutations(gene_seq)

        num_synonymous = sum(1 for m in mutations if m.mutation_type == 'synonymous')
        num_nonsynonymous = len(mutations) - num_synonymous

        logger.info(
            f"Variant analysis complete: {len(mutations)} mutations "
            f"({num_synonymous} synonymous, {num_nonsynonymous} non-synonymous)"
        )

        return {
            'gene_sequence': gene_seq,
            'protein_sequence': protein_seq,
            'mutations': mutations,
            'num_mutations': len(mutations),
            'num_synonymous': num_synonymous,
            'num_nonsynonymous': num_nonsynonymous,
        }
