import sqlite3
from typing import List, Tuple, Optional, Dict, Any, Iterable
from dataclasses import dataclass
import logging

from Bio.Seq import Seq
from Bio import Align

# ------------------------------------------------------------
# Logging (unchanged pattern)
# ------------------------------------------------------------
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# ------------------------------------------------------------
# Data classes
# ------------------------------------------------------------
@dataclass(frozen=True)
class MutationInfo:
    """
    Codon-level mutation call anchored to WT protein positions.

    position: 1-based WT amino-acid (codon) position.
    mutation_type: 'synonymous' or 'missense' (configurable: can emit legacy 'non-synonymous').
    """
    position: int
    wt_codon: str
    mut_codon: str
    wt_aa: str
    mut_aa: str
    mutation_type: str  # 'synonymous' or 'missense' (or legacy 'non-synonymous')

    def __str__(self) -> str:
        return f"{self.wt_aa}{self.position}{self.mut_aa} ({self.mutation_type})"


@dataclass(frozen=True)
class WTValidationMetadata:
    """
    Precomputed WT metadata required for deterministic extraction and downstream mapping.

    Required keys (assumed exact names):
      - start: int (0-based on plasmid string as provided)
      - strand: '+' or '-'
      - frame: 0/1/2
      - gene_length: int (nt)
      - match_type: str
      - identity: float 0..1
      - coverage: float 0..1
    """
    start: int
    strand: str
    frame: int
    gene_length: int
    match_type: str
    identity: float
    coverage: float

    @staticmethod
    def from_dict(d: Dict[str, Any]) -> "WTValidationMetadata":
        required = ("start", "strand", "frame", "gene_length", "match_type", "identity", "coverage")
        missing = [k for k in required if k not in d]
        if missing:
            raise KeyError(f"WT metadata missing required keys: {missing}")

        meta = WTValidationMetadata(
            start=int(d["start"]),
            strand=str(d["strand"]),
            frame=int(d["frame"]),
            gene_length=int(d["gene_length"]),
            match_type=str(d["match_type"]),
            identity=float(d["identity"]),
            coverage=float(d["coverage"]),
        )
        meta._validate()
        return meta

    def _validate(self) -> None:
        if self.strand not in ("+", "-"):
            raise ValueError(f"Invalid strand '{self.strand}'. Expected '+' or '-'.")
        if self.frame not in (0, 1, 2):
            raise ValueError(f"Invalid frame '{self.frame}'. Expected 0, 1, or 2.")
        if self.gene_length <= 0:
            raise ValueError("gene_length must be > 0.")


# ------------------------------------------------------------
# SequenceAnalyzer (modified)
# ------------------------------------------------------------
class SequenceAnalyzer:
    """
    Metadata-driven extraction + Biopython translation + Biopython protein alignment mapping.

    Key changes vs the original pipeline:
      - No WT discovery (_locate_wt_gene / Smith-Waterman) remains.
      - WT initialisation is driven by WTValidationMetadata and either WT plasmid or WT gene.
      - Mutation calling uses PairwiseAligner mapping between WT protein and translated variant protein.

    Notes on design:
      - Codon-level calls require the WT gene sequence (to know WT codons).
      - WT protein is used as the amino-acid anchor for mapping and reporting.
    """

    def __init__(
        self,
        wt_protein_seq: str,
        wt_metadata: Dict[str, Any],
        *,
        wt_plasmid_seq: Optional[str] = None,
        wt_gene_seq: Optional[str] = None,
        translation_table: Any = "Standard",
        aligner_mode: str = "global",           # 'global' (positional) or 'local' (robust to truncation/partials)
        match_score: float = 2.0,
        mismatch_score: float = -1.0,
        gap_score: float = -2.0,                # linear penalty via aligner.gap_score
        apply_frame_offset: bool = True,        # uses cds_start = start + frame if True
        strict: bool = True,                    # controls error vs warning for some issues
        mutation_label_style: str = "modern",   # 'modern' => missense; 'legacy' => non-synonymous
    ):
        self.strict = strict
        self.translation_table = translation_table
        self.apply_frame_offset = apply_frame_offset
        self.mutation_label_style = mutation_label_style

        self.wt_protein_seq = wt_protein_seq.strip().upper()
        self.wt_meta = WTValidationMetadata.from_dict(wt_metadata)

        # Store metadata fields with names similar to the original class, to reduce downstream breakage.
        self.wt_gene_start = self.wt_meta.start
        self.wt_gene_strand = self.wt_meta.strand
        self.wt_gene_frame = self.wt_meta.frame
        self.wt_gene_len = self.wt_meta.gene_length

        if wt_gene_seq is not None and wt_plasmid_seq is not None:
            raise ValueError("Provide either wt_gene_seq or wt_plasmid_seq, not both.")

        if wt_gene_seq is not None:
            self.wt_gene_seq = wt_gene_seq.strip().upper()
        elif wt_plasmid_seq is not None:
            self.wt_gene_seq = self._extract_gene_by_metadata(wt_plasmid_seq.strip().upper())
        else:
            raise ValueError("You must provide either wt_gene_seq or wt_plasmid_seq.")

        # Biopython PairwiseAligner setup (retained).
        self._aligner = Align.PairwiseAligner()
        if aligner_mode not in ("global", "local"):
            raise ValueError("aligner_mode must be 'global' or 'local'.")
        self._aligner.mode = aligner_mode
        self._aligner.match_score = float(match_score)
        self._aligner.mismatch_score = float(mismatch_score)
        self._aligner.gap_score = float(gap_score)

        self._sanity_check_wt_inputs()

        logger.info(
            "Initialised SequenceAnalyzer from WT metadata: start=%s strand=%s frame=%s gene_length=%s "
            "match_type=%s identity=%.3f coverage=%.3f aligner_mode=%s",
            self.wt_meta.start, self.wt_meta.strand, self.wt_meta.frame, self.wt_meta.gene_length,
            self.wt_meta.match_type, self.wt_meta.identity, self.wt_meta.coverage, aligner_mode,
        )

    # --- Alternative constructors (two requested patterns) ---
    @classmethod
    def from_wt_plasmid_and_metadata(
        cls,
        wt_plasmid_seq: str,
        wt_metadata: Dict[str, Any],
        wt_protein_seq: str,
        **kwargs: Any,
    ) -> "SequenceAnalyzer":
        return cls(
            wt_protein_seq=wt_protein_seq,
            wt_metadata=wt_metadata,
            wt_plasmid_seq=wt_plasmid_seq,
            wt_gene_seq=None,
            **kwargs,
        )

    @classmethod
    def from_wt_gene_and_metadata(
        cls,
        wt_gene_seq: str,
        wt_metadata: Dict[str, Any],
        wt_protein_seq: str,
        **kwargs: Any,
    ) -> "SequenceAnalyzer":
        return cls(
            wt_protein_seq=wt_protein_seq,
            wt_metadata=wt_metadata,
            wt_plasmid_seq=None,
            wt_gene_seq=wt_gene_seq,
            **kwargs,
        )

    # --- Internal helpers: extraction, translation, alignment mapping ---
    def _cds_start(self, plasmid_len: int) -> int:
        """
        Compute coding-start coordinate in plasmid coordinates.

        Assumed semantics (explicit):
          cds_start = (start + frame) mod plasmid_len   if apply_frame_offset
          cds_start = start mod plasmid_len             otherwise
        """
        if plasmid_len <= 0:
            raise ValueError("Plasmid length must be > 0.")
        s = self.wt_meta.start % plasmid_len
        if self.apply_frame_offset:
            s = (s + self.wt_meta.frame) % plasmid_len
        return s

    def _extract_circular(self, seq: str, start: int, length: int) -> str:
        """
        Extract length bases from circular sequence starting at start (0-based).
        """
        n = len(seq)
        if n == 0:
            raise ValueError("Sequence is empty.")
        if length <= 0:
            raise ValueError("length must be > 0.")
        if length > n:
            raise ValueError(f"Cannot extract length={length} from sequence length={n}.")
        start = start % n
        end = start + length
        if end <= n:
            return seq[start:end]
        return seq[start:] + seq[: (end % n)]

    def _extract_gene_by_metadata(self, plasmid_seq: str) -> str:
        """
        Extraction using metadata: start/frame/gene_length + strand.
        Returns gene in CODING ORIENTATION.
        """
        plasmid_len = len(plasmid_seq)
        cds_start = self._cds_start(plasmid_len)
        window = self._extract_circular(plasmid_seq, cds_start, self.wt_gene_len)

        if self.wt_gene_strand == "-":
            window = str(Seq(window).reverse_complement())
        return window

    def _translate_dna(self, dna_seq: str) -> str:
        """
        Translate DNA to protein using Biopython Seq.translate.

        Implementation uses to_stop=False to preserve 1 codon => 1 residue mapping
        (stop codons appear as '*'), which simplifies codon-index mapping.
        """
        dna_seq = dna_seq.strip().upper()

        # Enforce/handle length multiple of 3 for stable codon indexing.
        if len(dna_seq) % 3 != 0:
            msg = f"DNA length {len(dna_seq)} is not divisible by 3."
            if self.strict:
                raise ValueError(msg)
            else:
                logger.warning(msg + " Trimming to full codons for translation.")
                dna_seq = dna_seq[: (len(dna_seq) // 3) * 3]

        try:
            return str(Seq(dna_seq).translate(table=self.translation_table, to_stop=False, cds=False))
        except Exception as e:
            # Biopython raises for invalid codons; in strict mode fail fast.
            if self.strict:
                raise
            logger.warning("Translation failed (non-strict). Returning empty protein. Error=%r", e)
            return ""

    def _normalise_protein_for_alignment(self, protein: str) -> str:
        """
        Normalise proteins for alignment.
        - Upper-case
        - Remove trailing '*' since WT protein inputs typically omit terminal stop symbols.
        """
        p = protein.strip().upper()
        return p.rstrip("*")

    def _align_wt_to_variant(self, wt_protein: str, variant_protein: str):
        """
        Align WT protein (target) to variant protein (query) using PairwiseAligner.
        """
        if not wt_protein or not variant_protein:
            return None
        alignments = self._aligner.align(wt_protein, variant_protein)
        if len(alignments) == 0:
            return None
        return alignments[0]

    def _iter_aligned_index_pairs(self, alignment) -> Iterable[Tuple[int, int]]:
        """
        Yield (wt_res_idx, var_res_idx) for residue positions that are aligned to each other.

        Uses alignment.aligned blocks:
          alignment.aligned[0] => blocks on target (WT)
          alignment.aligned[1] => blocks on query (variant)
        """
        # alignment.aligned is a NumPy array in Biopython; handle generically.
        blocks_wt = alignment.aligned[0]
        blocks_var = alignment.aligned[1]

        # Convert to Python lists of pairs (optional safety)
        try:
            blocks_wt = blocks_wt.tolist()
            blocks_var = blocks_var.tolist()
        except Exception:
            pass

        for (w_s, w_e), (v_s, v_e) in zip(blocks_wt, blocks_var):
            length = min(int(w_e) - int(w_s), int(v_e) - int(v_s))
            for k in range(length):
                yield int(w_s) + k, int(v_s) + k

    def _translate_codon(self, codon: str) -> str:
        """
        Translate a single codon with Biopython, returning a 1-letter amino acid or 'X' on failure.
        """
        codon = codon.upper()
        if len(codon) != 3:
            return "X"
        try:
            return str(Seq(codon).translate(table=self.translation_table, to_stop=False, cds=False))
        except Exception:
            return "X"

    def _sanity_check_wt_inputs(self) -> None:
        """
        Lightweight checks for common metadata/input mismatches.
        """
        if len(self.wt_gene_seq) != self.wt_gene_len:
            msg = (
                f"WT gene length mismatch: len(wt_gene_seq)={len(self.wt_gene_seq)} vs "
                f"metadata gene_length={self.wt_gene_len}"
            )
            if self.strict:
                raise ValueError(msg)
            logger.warning(msg)

        # Optional: Confirm WT gene translation aligns reasonably to WT protein.
        wt_gene_protein = self._normalise_protein_for_alignment(self._translate_dna(self.wt_gene_seq))
        wt_protein_norm = self._normalise_protein_for_alignment(self.wt_protein_seq)

        alignment = self._align_wt_to_variant(wt_protein_norm, wt_gene_protein)
        if alignment is None:
            msg = "Could not align WT gene translation to WT protein; check inputs or scoring."
            if self.strict:
                raise ValueError(msg)
            logger.warning(msg)
            return

        # crude identity/coverage (for diagnostics only; metadata identity/coverage may come from upstream)
        aligned_len = 0
        matches = 0
        for w_i, g_i in self._iter_aligned_index_pairs(alignment):
            aligned_len += 1
            if wt_protein_norm[w_i] == wt_gene_protein[g_i]:
                matches += 1
        identity = (matches / aligned_len) if aligned_len else 0.0
        coverage = (aligned_len / len(wt_protein_norm)) if wt_protein_norm else 0.0

        # Compare to metadata heuristically; do not hard fail unless strict and wildly off.
        if self.strict and (coverage + 1e-9) < min(0.8, self.wt_meta.coverage):
            raise ValueError(
                f"WT protein coverage from alignment appears low (observed≈{coverage:.3f}, "
                f"metadata={self.wt_meta.coverage:.3f}). Check metadata coordinate semantics."
            )
        if (identity + 1e-9) < min(0.8, self.wt_meta.identity):
            logger.warning(
                "WT protein identity appears lower than metadata (observed≈%.3f, metadata=%.3f).",
                identity, self.wt_meta.identity,
            )

    # ------------------------------------------------------------
    # Public API (variant processing)
    # ------------------------------------------------------------
    def extract_gene_from_variant(self, variant_plasmid_seq: str) -> str:
        """
        Extract gene from variant plasmid using WT metadata.

        Required behaviour:
          - uses metadata start, frame, strand, gene_length;
          - handles circular plasmids;
          - returns gene in coding orientation (reverse-complement if strand == '-').

        Note:
          - This assumes the variant plasmid matches WT coordinate system (no upstream large indels).
        """
        variant_plasmid_seq = variant_plasmid_seq.strip().upper()
        return self._extract_gene_by_metadata(variant_plasmid_seq)

    def identify_mutations(self, variant_gene_seq: str) -> List[MutationInfo]:
        """
        Call synonymous vs missense mutations by:
          1) translating variant DNA (Biopython);
          2) aligning translated variant protein to WT protein (Biopython PairwiseAligner);
          3) using alignment mapping to compare WT vs variant codons at aligned residue indices.

        Only substitutions at aligned residue positions are called.
        Insertions/deletions (gaps) are not classified into synonymous/missense by this function.
        """
        variant_gene_seq = variant_gene_seq.strip().upper()

        # Translate variant gene to protein (full-length, includes '*' if internal stops)
        variant_protein_full = self._translate_dna(variant_gene_seq)
        wt_protein_norm = self._normalise_protein_for_alignment(self.wt_protein_seq)
        var_protein_norm = self._normalise_protein_for_alignment(variant_protein_full)

        alignment = self._align_wt_to_variant(wt_protein_norm, var_protein_norm)
        if alignment is None:
            logger.warning("No protein alignment found; returning no mutation calls.")
            return []

        # Mutation calling
        muts: List[MutationInfo] = []
        for wt_idx, var_idx in self._iter_aligned_index_pairs(alignment):
            # Map residue indices to codon indices (1 residue per codon, by design)
            wt_codon_off = wt_idx * 3
            var_codon_off = var_idx * 3
            if wt_codon_off + 3 > len(self.wt_gene_seq) or var_codon_off + 3 > len(variant_gene_seq):
                continue

            wt_codon = self.wt_gene_seq[wt_codon_off:wt_codon_off + 3]
            mut_codon = variant_gene_seq[var_codon_off:var_codon_off + 3]

            if wt_codon == mut_codon:
                continue

            # Determine amino acids:
            # - WT AA from WT protein (alignment target)
            # - Variant AA from translated variant protein (alignment query)
            wt_aa = wt_protein_norm[wt_idx] if wt_idx < len(wt_protein_norm) else self._translate_codon(wt_codon)
            mut_aa = var_protein_norm[var_idx] if var_idx < len(var_protein_norm) else self._translate_codon(mut_codon)

            if wt_aa == mut_aa:
                label = "synonymous"
            else:
                label = "missense" if self.mutation_label_style == "modern" else "non-synonymous"

            muts.append(
                MutationInfo(
                    position=wt_idx + 1,
                    wt_codon=wt_codon,
                    mut_codon=mut_codon,
                    wt_aa=wt_aa,
                    mut_aa=mut_aa,
                    mutation_type=label,
                )
            )

        return muts

    def analyze_variant(self, variant_plasmid_seq: str) -> dict:
        """
        Full analysis returning dict compatible with existing SequenceDatabase storage:
          - gene_sequence
          - protein_sequence
          - mutations (List[MutationInfo])
          - num_mutations, num_synonymous, num_nonsynonymous
        """
        gene_seq = self.extract_gene_from_variant(variant_plasmid_seq)
        protein_seq = self._normalise_protein_for_alignment(self._translate_dna(gene_seq))
        mutations = self.identify_mutations(gene_seq)

        num_syn = sum(1 for m in mutations if m.mutation_type == "synonymous")
        num_non = len(mutations) - num_syn

        logger.info("Variant analysis complete: %s mutations (%s syn, %s non-syn/missense)", len(mutations), num_syn, num_non)

        return {
            "gene_sequence": gene_seq,
            "protein_sequence": protein_seq,
            "mutations": mutations,
            "num_mutations": len(mutations),
            "num_synonymous": num_syn,
            "num_nonsynonymous": num_non,
        }


# ------------------------------------------------------------
# SequenceDatabase (preserved logic; unchanged schema)
# ------------------------------------------------------------
class SequenceDatabase:
    """Database handler for storing sequence analysis results."""

    def __init__(self, db_path: str = "directed_evolution.db"):
        self.db_path = db_path
        self.conn = sqlite3.connect(self.db_path)
        self.cursor = self.conn.cursor()
        self._create_tables()
        logger.info(f"Connected to database: {self.db_path}")

    def _create_tables(self):
        self.cursor.execute("""
            CREATE TABLE IF NOT EXISTS variants (
                variant_id INTEGER PRIMARY KEY,
                gene_sequence TEXT,
                protein_sequence TEXT,
                num_mutations INTEGER,
                num_synonymous INTEGER,
                num_nonsynonymous INTEGER
            )
        """)

        self.cursor.execute("""
            CREATE TABLE IF NOT EXISTS mutations (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                variant_id INTEGER,
                position INTEGER,
                wt_codon TEXT,
                mut_codon TEXT,
                wt_aa TEXT,
                mut_aa TEXT,
                mutation_type TEXT,
                FOREIGN KEY (variant_id) REFERENCES variants(variant_id)
            )
        """)

        self.conn.commit()

    def store_variant(self, variant_id: int, analysis: dict):
        self.cursor.execute("""
            INSERT OR REPLACE INTO variants
            (variant_id, gene_sequence, protein_sequence, num_mutations,
             num_synonymous, num_nonsynonymous)
            VALUES (?, ?, ?, ?, ?, ?)
        """, (
            variant_id,
            analysis["gene_sequence"],
            analysis["protein_sequence"],
            analysis["num_mutations"],
            analysis["num_synonymous"],
            analysis["num_nonsynonymous"],
        ))

        for mutation in analysis["mutations"]:
            self.cursor.execute("""
                INSERT INTO mutations
                (variant_id, position, wt_codon, mut_codon, wt_aa, mut_aa, mutation_type)
                VALUES (?, ?, ?, ?, ?, ?, ?)
            """, (
                variant_id,
                mutation.position,
                mutation.wt_codon,
                mutation.mut_codon,
                mutation.wt_aa,
                mutation.mut_aa,
                mutation.mutation_type,
            ))

        self.conn.commit()
        logger.info(f"Stored analysis for variant {variant_id}")

    def get_variant(self, variant_id: int) -> Optional[dict]:
        self.cursor.execute("SELECT * FROM variants WHERE variant_id = ?", (variant_id,))
        row = self.cursor.fetchone()
        if not row:
            return None

        self.cursor.execute("SELECT * FROM mutations WHERE variant_id = ?", (variant_id,))
        mutations = self.cursor.fetchall()

        return {
            "variant_id": row[0],
            "gene_sequence": row[1],
            "protein_sequence": row[2],
            "num_mutations": row[3],
            "num_synonymous": row[4],
            "num_nonsynonymous": row[5],
            "mutations": mutations,
        }

    def close(self):
        self.conn.close()
        logger.info("Database connection closed")


# ------------------------------------------------------------
# Example instantiation patterns (requested) + unit-test style demos
# ------------------------------------------------------------
def _example_metadata_forward() -> Dict[str, Any]:
    return {
        "start": 6,                 # anchor in plasmid coordinates
        "strand": "+",
        "frame": 0,                 # cds_start = start + frame
        "gene_length": 9,           # 3 codons
        "match_type": "validated",
        "identity": 1.0,
        "coverage": 1.0,
    }

def demo_instantiation_patterns_and_calls():
    # WT gene encodes MKF: ATG AAA TTT
    wt_gene = "ATGAAATTT"
    wt_protein = "MKF"

    # Pattern (a): WT plasmid + metadata + WT protein
    wt_plasmid = "AAACCC" + wt_gene + "GGGTTT"  # gene begins at index 6, no wrap
    analyzer_a = SequenceAnalyzer.from_wt_plasmid_and_metadata(
        wt_plasmid_seq=wt_plasmid,
        wt_metadata=_example_metadata_forward(),
        wt_protein_seq=wt_protein,
        strict=True,
        aligner_mode="global",
        mutation_label_style="modern",  # emits 'missense'
    )

    # Pattern (b): Pre-extracted WT gene + metadata + WT protein
    analyzer_b = SequenceAnalyzer.from_wt_gene_and_metadata(
        wt_gene_seq=wt_gene,
        wt_metadata=_example_metadata_forward(),
        wt_protein_seq=wt_protein,
        strict=True,
        aligner_mode="global",
        mutation_label_style="modern",
    )

    # Variant 1: AAA -> AAG (synonymous K)
    var1_plasmid = "AAACCC" + "ATGAAGTTT" + "GGGTTT"
    res1 = analyzer_a.analyze_variant(var1_plasmid)
    assert res1["num_mutations"] == 1
    assert res1["mutations"][0].mutation_type == "synonymous"
    assert str(res1["mutations"][0]) == "K2K (synonymous)"

    # Variant 2: TTT -> TAT (missense F->Y)
    var2_plasmid = "AAACCC" + "ATGAAATAT" + "GGGTTT"
    res2 = analyzer_b.analyze_variant(var2_plasmid)
    assert res2["num_mutations"] == 1
    assert res2["mutations"][0].mutation_type == "missense"
    assert str(res2["mutations"][0]) == "F3Y (missense)"

def test_reverse_strand_extraction_and_synonymous():
    # WT gene: ATG AAA TTT => MKF
    wt_gene = "ATGAAATTT"
    wt_protein = "MKF"

    # Reverse complement of WT gene is AAATTTCAT
    rc_gene = str(Seq(wt_gene).reverse_complement())

    # Put reverse-strand gene window into plasmid at index 3
    wt_plasmid = "GGG" + rc_gene + "CCC"
    meta = {
        "start": 3,
        "strand": "-",
        "frame": 0,
        "gene_length": 9,
        "match_type": "validated",
        "identity": 1.0,
        "coverage": 1.0,
    }

    analyzer = SequenceAnalyzer.from_wt_plasmid_and_metadata(
        wt_plasmid_seq=wt_plasmid,
        wt_metadata=meta,
        wt_protein_seq=wt_protein,
        strict=True,
        aligner_mode="global",
        mutation_label_style="modern",
    )

    # Confirm WT gene recovered in coding orientation
    assert analyzer.wt_gene_seq == wt_gene

    # Make a synonymous change in the third codon: TTT -> TTC (F stays F)
    var_gene = "ATGAAATTC"
    var_plasmid = "GGG" + str(Seq(var_gene).reverse_complement()) + "CCC"

    res = analyzer.analyze_variant(var_plasmid)
    assert res["num_mutations"] == 1
    assert res["mutations"][0].mutation_type == "synonymous"
    assert str(res["mutations"][0]) == "F3F (synonymous)"
