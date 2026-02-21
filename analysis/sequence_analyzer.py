import logging
from typing import List, Optional
from Bio import Align
from Bio.Seq import Seq
from .mutation_models import MutationInfo

logger = logging.getLogger(__name__)


class SequenceAnalyzer:

    def __init__(self, wt_protein_seq: str, wt_plasmid_seq: str):
        self.wt_protein_seq = wt_protein_seq.upper().strip()
        self.wt_plasmid_seq = wt_plasmid_seq.upper().strip()

        self.wt_gene_seq: Optional[str] = None
        self.wt_gene_start: Optional[int] = None
        self.wt_gene_strand: Optional[str] = None

        self._aligner = Align.PairwiseAligner()
        self._aligner.mode = "local"
        self._aligner.match_score = 2.0
        self._aligner.mismatch_score = -1.0
        self._aligner.gap_score = -1.0

        self._locate_wt_gene()

    def _translate(self, dna_seq: str) -> str:
        return str(Seq(dna_seq).translate(to_stop=True))

    def _reverse_complement(self, seq: str) -> str:
        return str(Seq(seq).reverse_complement())

    def _smith_waterman(self, seq1, seq2):
        alignments = self._aligner.align(seq1, seq2)
        if not alignments:
            return 0.0, 0, 0

        alignment = alignments[0]
        blocks = alignment.aligned[0]
        start = blocks[0][0]
        end = blocks[-1][1]
        return alignment.score, start, end

    def _locate_wt_gene(self):
        plasmid_len = len(self.wt_plasmid_seq)
        extended = self.wt_plasmid_seq + self.wt_plasmid_seq

        best_score = float("-inf")
        best = None

        for strand, seq in [("+", extended),
                            ("-", self._reverse_complement(extended))]:

            for frame in range(3):
                translated = self._translate(seq[frame:])
                score, aa_start, aa_end = self._smith_waterman(
                    translated, self.wt_protein_seq
                )

                if score > best_score:
                    dna_start = frame + aa_start * 3
                    dna_end = frame + aa_end * 3
                    gene_seq = seq[dna_start:dna_end]

                    best_score = score
                    best = (strand, gene_seq, dna_start % plasmid_len)

        if not best or best_score < len(self.wt_protein_seq) * 1.5:
            raise ValueError("Could not locate WT gene in plasmid")

        strand, gene_seq, start = best

        # Validation
        if self._translate(gene_seq) != self.wt_protein_seq:
            raise ValueError("WT gene translation does not match WT protein")

        self.wt_gene_seq = gene_seq
        self.wt_gene_start = start
        self.wt_gene_strand = strand

    def extract_gene_from_variant(self, variant_plasmid_seq: str) -> str:
        variant_plasmid_seq = variant_plasmid_seq.upper().strip()
        gene_len = len(self.wt_gene_seq)

        if self.wt_gene_start + gene_len <= len(variant_plasmid_seq):
            segment = variant_plasmid_seq[
                self.wt_gene_start:self.wt_gene_start + gene_len
            ]
        else:
            segment = (
                variant_plasmid_seq[self.wt_gene_start:] +
                variant_plasmid_seq[:gene_len - (
                    len(variant_plasmid_seq) - self.wt_gene_start
                )]
            )

        if self.wt_gene_strand == "-":
            segment = self._reverse_complement(segment)

        return segment

    def identify_mutations(self, variant_gene: str) -> List[MutationInfo]:
        mutations = []

        for i in range(0, len(self.wt_gene_seq), 3):
            wt_codon = self.wt_gene_seq[i:i+3]
            mut_codon = variant_gene[i:i+3]

            if wt_codon != mut_codon:
                wt_aa = self._translate(wt_codon)
                mut_aa = self._translate(mut_codon)

                mutation_type = (
                    "synonymous" if wt_aa == mut_aa else "non-synonymous"
                )

                mutations.append(
                    MutationInfo(
                        position=(i // 3) + 1,
                        wt_codon=wt_codon,
                        mut_codon=mut_codon,
                        wt_aa=wt_aa,
                        mut_aa=mut_aa,
                        mutation_type=mutation_type,
                    )
                )

        return mutations

    def analyze_variant(self, variant_plasmid: str):
        gene_seq = self.extract_gene_from_variant(variant_plasmid)
        protein_seq = self._translate(gene_seq)
        mutations = self.identify_mutations(gene_seq)

        return {
            "gene_sequence": gene_seq,
            "protein_sequence": protein_seq,
            "mutations": mutations,
            "num_mutations": len(mutations),
            "num_synonymous": sum(m.mutation_type == "synonymous" for m in mutations),
            "num_nonsynonymous": sum(m.mutation_type == "non-synonymous" for m in mutations),
        }