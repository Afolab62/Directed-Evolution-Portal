import logging
from typing import List, Tuple, Optional

from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

try:
    from database.models import MutationInfo, StagingResult
except ImportError:  # pragma: no cover
    from models import MutationInfo, StagingResult


logger = logging.getLogger(__name__)


class SequenceAnalyzer:
    """
    Sequence analysis and mutation detection.

    Gene coordinates are provided by a StagingResult from the plasmid
    validation pipeline rather than being re-derived at runtime, so
    _locate_wt_gene() is only called as a fallback when no staging result
    is available.
    """

    # Standard genetic code (used for per-codon lookup in mutation identification)
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
            When provided (and is_valid=True) this is used directly, skipping
            the Smith-Waterman search entirely. Pass None to fall back to the
            built-in alignment.
        """
        self.wt_protein_seq = wt_protein_seq.upper().strip()
        self.wt_plasmid_seq = wt_plasmid_seq.upper().strip()
        self.wt_gene_seq: Optional[str] = None
        self.wt_gene_start: Optional[int] = None
        self.staging_result = staging_result

        if staging_result is not None:
            self._load_from_staging(staging_result)
        else:
            self._locate_wt_gene()

    # ------------------------------------------------------------------
    # Staging result integration
    # ------------------------------------------------------------------

    def _load_from_staging(self, sr: StagingResult) -> None:
        """
        Populate wt_gene_seq / wt_gene_start directly from a StagingResult
        row produced by the plasmid validation pipeline.
        """
        if not sr.is_valid:
            raise ValueError(
                f"staging_result id={sr.id} is marked is_valid=False "
                f"(notes: {sr.notes}). Cannot proceed with this result."
            )

        plasmid = self.wt_plasmid_seq
        start = sr.start_nt
        end = sr.end_nt_exclusive

        # Extract gene, honouring origin-spanning genes
        if sr.wraps_origin:
            gene_seq = plasmid[start:] + plasmid[:end % len(plasmid)]
        else:
            gene_seq = plasmid[start:end]

        # Reverse-complement if the gene sits on the minus strand
        if sr.strand == '-':
            gene_seq = self._reverse_complement(gene_seq)

        # Shift into the correct reading frame then trim to codon boundary
        gene_seq = gene_seq[sr.frame:]
        gene_seq = gene_seq[:len(gene_seq) - len(gene_seq) % 3]

        self.wt_gene_seq = gene_seq
        self.wt_gene_start = start

        logger.info(
            f"Gene loaded from staging_result id={sr.id}: "
            f"strand={sr.strand}, frame={sr.frame}, start={start}, "
            f"length={len(gene_seq)}, identity={sr.identity:.3f}, "
            f"coverage={sr.coverage:.3f}"
        )

    # ------------------------------------------------------------------
    # Internal alignment helpers (fallback only)
    # ------------------------------------------------------------------

    def _smith_waterman(self, seq1: str, seq2: str, match=2, mismatch=-1, gap=-1) -> Tuple[int, int, int]:
        """
        Smith-Waterman local alignment using BioPython's PairwiseAligner.
        Returns: (score, start_pos_in_seq1, end_pos_in_seq1)
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
        """Translate DNA to protein sequence using BioPython."""
        trimmed = dna_seq[:len(dna_seq) - len(dna_seq) % 3]
        return str(Seq(trimmed).translate(to_stop=True))

    def _reverse_complement(self, seq: str) -> str:
        """Generate reverse complement of a DNA sequence using BioPython."""
        return str(Seq(seq).reverse_complement())

    def _locate_wt_gene(self) -> None:
        """
        Fallback: locate gene via Smith-Waterman on translated sequences.
        Called only when no StagingResult is supplied.
        """
        logger.info("No staging result provided — locating WT gene via Smith-Waterman alignment...")

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
            raise ValueError("Could not locate WT gene in plasmid!")

        self.wt_gene_seq = best_match['gene_seq']
        self.wt_gene_start = best_match['start']

        logger.info(
            f"WT gene found: score={best_match['score']}, "
            f"start={self.wt_gene_start}, length={len(self.wt_gene_seq)}"
        )

    # ------------------------------------------------------------------
    # Core analysis
    # ------------------------------------------------------------------

    def extract_gene_from_variant(self, variant_plasmid_seq: str) -> str:
        """Extract the gene region from a variant plasmid using WT coordinates."""
        variant_plasmid_seq = variant_plasmid_seq.upper().strip()
        gene_len = len(self.wt_gene_seq)

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

        # Apply strand / frame correction from staging if available
        if self.staging_result is not None:
            if self.staging_result.strand == '-':
                gene_seq = self._reverse_complement(gene_seq)
            gene_seq = gene_seq[self.staging_result.frame:]
            gene_seq = gene_seq[:len(gene_seq) - len(gene_seq) % 3]

        return gene_seq

    def identify_mutations(self, variant_gene: str) -> List[MutationInfo]:
        """Identify mutations by comparing variant gene to WT codon by codon."""
        mutations = []

        if len(variant_gene) != len(self.wt_gene_seq):
            raise ValueError(
                f"Gene length mismatch (variant={len(variant_gene)} bp, WT={len(self.wt_gene_seq)} bp). "
                "This likely indicates an indel or frameshift. Codon-level comparison cannot proceed."
            )

        for i in range(0, len(self.wt_gene_seq) - 2, 3):
            wt_codon = self.wt_gene_seq[i:i+3]
            mut_codon = variant_gene[i:i+3]

            if wt_codon != mut_codon:
                wt_aa = self.GENETIC_CODE.get(wt_codon, 'X')
                mut_aa = self.GENETIC_CODE.get(mut_codon, 'X')
                mutations.append(MutationInfo(
                    position=(i // 3) + 1,
                    wt_codon=wt_codon,
                    mut_codon=mut_codon,
                    wt_aa=wt_aa,
                    mut_aa=mut_aa,
                    mutation_type='synonymous' if wt_aa == mut_aa else 'non-synonymous'
                ))

        return mutations

    def analyze_variant(self, variant_plasmid: str) -> dict:
        """
        Complete analysis of a variant plasmid.
        Returns a dict with gene_seq, protein_seq, and mutation details.
        """
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
            'num_nonsynonymous': num_nonsynonymous,
        }
