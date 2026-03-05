"""
Service for sequence processing: DNA->Protein translation and mutation detection.
Uses Smith-Waterman local alignment to robustly locate the gene of interest in
circular plasmid DNA, then locks those coordinates for all variant extractions.
"""

from typing import Dict, List, Tuple, Optional
import logging

logger = logging.getLogger(__name__)


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
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def _clean(seq: str) -> str:
    return seq.strip().upper().replace("\n", "").replace("\r", "")


def _translate(dna_seq: str) -> str:
    """Translate DNA to protein, stopping at the first stop codon."""
    protein = []
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        aa = GENETIC_CODE.get(codon, 'X')
        if aa == '*':
            break
        protein.append(aa)
    return ''.join(protein)


def _translate_full(dna_seq: str) -> str:
    """Translate DNA to protein including stop codons (as '*')."""
    protein = []
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        aa = GENETIC_CODE.get(codon, 'X')
        protein.append(aa)
    return ''.join(protein)


def _reverse_complement(seq: str) -> str:
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement.get(b, b) for b in reversed(seq))


def _needleman_wunsch(seq_a: str, seq_b: str,
                      match: int = 2, mismatch: int = -1, gap: int = -2) -> Tuple[str, str]:
    """
    Needleman-Wunsch global alignment of two protein sequences.
    Returns (aligned_seq_a, aligned_seq_b) with '-' gap characters.
    Used to compute alignment-aware positions for mutations near indels.
    """
    n, m = len(seq_a), len(seq_b)
    score = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[0] * (m + 1) for _ in range(n + 1)]  # 0=diag 1=up 2=left

    for i in range(1, n + 1):
        score[i][0] = i * gap
        trace[i][0] = 1
    for j in range(1, m + 1):
        score[0][j] = j * gap
        trace[0][j] = 2

    for i in range(1, n + 1):
        a = seq_a[i - 1]
        for j in range(1, m + 1):
            b = seq_b[j - 1]
            s_diag = score[i-1][j-1] + (match if a == b else mismatch)
            s_up   = score[i-1][j] + gap
            s_left = score[i][j-1] + gap
            best = s_diag; t = 0
            if s_up   > best: best = s_up;   t = 1
            if s_left > best: best = s_left; t = 2
            score[i][j] = best
            trace[i][j] = t

    aligned_a: List[str] = []
    aligned_b: List[str] = []
    i, j = n, m
    while i > 0 or j > 0:
        t = trace[i][j] if i >= 0 and j >= 0 else 0
        if i > 0 and j > 0 and t == 0:
            aligned_a.append(seq_a[i-1]); aligned_b.append(seq_b[j-1]); i -= 1; j -= 1
        elif i > 0 and (j == 0 or t == 1):
            aligned_a.append(seq_a[i-1]); aligned_b.append('-'); i -= 1
        else:
            aligned_a.append('-'); aligned_b.append(seq_b[j-1]); j -= 1
    aligned_a.reverse(); aligned_b.reverse()
    return ''.join(aligned_a), ''.join(aligned_b)


def _build_wt_position_map(aligned_wt: str, aligned_var: str) -> Dict[int, int]:
    """Map each 1-based WT residue position to its alignment column index."""
    wt_map: Dict[int, int] = {}
    wt_pos = 0
    for aln_pos, (wa, _) in enumerate(zip(aligned_wt, aligned_var), start=1):
        if wa != '-':
            wt_pos += 1
            wt_map[wt_pos] = aln_pos
    return wt_map


def _estimate_rotation_offset(wt_plasmid: str, variant_plasmid: str) -> Optional[int]:
    """
    Estimate the circular rotation offset (in nt) between WT and a variant
    plasmid assembly.  High-throughput sequencing assemblers may produce
    assemblies that start at a different position on the circle, shifting the
    gene coordinates by a fixed amount.

    Short anchors from the WT are located in the variant; the most-voted
    offset is returned.  Returns None if no consistent offset is found.
    """
    wt  = _clean(wt_plasmid)
    var = _clean(variant_plasmid)
    length = len(wt)
    anchor_positions = list(range(0, length, 350))

    for k in (150, 120, 100, 80, 60, 40):
        votes: Dict[int, int] = {}
        for pos in anchor_positions:
            if pos + k > length:
                continue
            anchor = wt[pos: pos + k]
            idx = var.find(anchor)
            while idx >= 0:
                offset = (pos - idx) % length
                votes[offset] = votes.get(offset, 0) + 1
                idx = var.find(anchor, idx + 1)
        if votes:
            return max(votes.items(), key=lambda kv: (kv[1], -kv[0]))[0]
    return None


def _smith_waterman(seq1: str, seq2: str, match=2, mismatch=-1, gap=-1) -> Tuple[int, int, int]:
    """
    Smith-Waterman local alignment.
    Returns (score, start_in_seq1, end_in_seq1).
    """
    m, n = len(seq1), len(seq2)
    H = [[0] * (n + 1) for _ in range(m + 1)]
    max_score = 0
    max_pos = (0, 0)

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diag = H[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            delete = H[i-1][j] + gap
            insert = H[i][j-1] + gap
            H[i][j] = max(0, diag, delete, insert)
            if H[i][j] > max_score:
                max_score = H[i][j]
                max_pos = (i, j)

    # Traceback to find start
    i, j = max_pos
    while i > 0 and j > 0 and H[i][j] > 0:
        i -= 1
        j -= 1

    return max_score, i, max_pos[0]


def locate_gene_sw(wt_plasmid_seq: str, wt_protein_seq: str) -> Tuple[str, int, int, bool]:
    """
    Use Smith-Waterman alignment to locate the gene encoding wt_protein_seq
    in the plasmid (handles circular DNA by doubling the sequence).

    Returns (gene_dna, gene_start, gene_length, is_rc_strand).
    gene_dna is in 5'->3' coding orientation.
    gene_start / gene_length are coordinates on the original plus-strand plasmid.
    """
    plasmid = _clean(wt_plasmid_seq)
    protein = _clean(wt_protein_seq)
    L = len(plasmid)

    # Double plasmid to handle genes that wrap around the origin
    extended = plasmid + plasmid

    best_score = 0
    best_gene_seq = None
    best_start = 0
    best_length = 0
    best_is_rc = False

    for strand_seq, is_rc in [(extended, False), (_reverse_complement(extended), True)]:
        for frame in range(3):
            translated = _translate_full(strand_seq[frame:])
            score, start_aa, end_aa = _smith_waterman(translated, protein)

            if score > best_score:
                dna_start = frame + start_aa * 3
                dna_end = frame + end_aa * 3
                gene_dna = strand_seq[dna_start:dna_end]  # coding-orientation DNA

                best_score = score
                # When is_rc=True, dna_start is a position in rc(extended).
                # The equivalent start on the plus-strand extended is (2L - dna_end).
                if is_rc:
                    best_start = (2 * L - dna_end) % L
                else:
                    best_start = dna_start % L
                best_length = dna_end - dna_start
                best_gene_seq = gene_dna
                best_is_rc = is_rc

    min_acceptable = len(protein) * 1.5  # need at least 75% identity
    if best_score < min_acceptable or best_gene_seq is None:
        raise ValueError(
            f"Could not reliably locate WT gene in plasmid "
            f"(best score {best_score:.0f}, threshold {min_acceptable:.0f})"
        )

    gene_protein = _translate(best_gene_seq)
    print(f"  Gene located: start={best_start}, length={best_length}bp ({len(gene_protein)} AA), score={best_score:.0f}")
    print(f"  First 50 AA: {gene_protein[:50]}")
    print(f"  Wraps origin: {best_start + best_length > L}")
    print(f"  Strand: {'minus (RC)' if best_is_rc else 'plus'}")

    return best_gene_seq, best_start, best_length, best_is_rc


def extract_gene(plasmid_seq: str, gene_start: int, gene_length: int) -> str:
    """Extract gene_length bp starting at gene_start from a circular plasmid."""
    seq = _clean(plasmid_seq)
    L = len(seq)
    if gene_start + gene_length <= L:
        return seq[gene_start:gene_start + gene_length]
    # Wraps around origin
    part1 = seq[gene_start:]
    part2 = seq[:gene_length - len(part1)]
    return part1 + part2


def identify_mutations(wt_gene_dna: str, var_gene_dna: str) -> List[Dict]:
    """
    Codon-by-codon comparison.  Returns list of mutation dicts with keys:
        position, wt_aa, mut_aa, wt_codon, mut_codon, mutation_type, aa_change

    Field names are standardised to match mutation_analysis.py so that
    experiments.py can consume results from either code path without branching.
    aligned_position is added downstream by analyze_variant_batch() after the
    global protein alignment step.
    """
    mutations = []

    if len(wt_gene_dna) != len(var_gene_dna):
        logger.warning(
            f"Gene length mismatch: WT={len(wt_gene_dna)}, variant={len(var_gene_dna)}. "
            "Truncating to shorter length."
        )

    n_codons = min(len(wt_gene_dna), len(var_gene_dna)) // 3

    for i in range(n_codons):
        s = i * 3
        wt_codon  = wt_gene_dna[s:s+3]
        var_codon = var_gene_dna[s:s+3]

        if wt_codon == var_codon:
            continue

        wt_aa  = GENETIC_CODE.get(wt_codon,  'X')
        var_aa = GENETIC_CODE.get(var_codon, 'X')

        # Skip stop codons mid-sequence
        if wt_aa == '*' or var_aa == '*':
            continue

        pos = i + 1  # 1-based
        mutations.append({
            'position':      pos,
            'wt_aa':         wt_aa,
            'mut_aa':        var_aa,
            'wt_codon':      wt_codon,
            'mut_codon':     var_codon,
            'mutation_type': 'synonymous' if wt_aa == var_aa else 'non-synonymous',
            'aa_change':     f'{wt_aa}{pos}{var_aa}',
        })

    return mutations



def locate_gene_fast(wt_plasmid_seq: str, wt_protein_seq: str) -> Tuple[str, int, int, bool]:
    """
    Quickly locate the gene encoding wt_protein_seq in the plasmid.

    Strategy:
    1. Translate every frame (0-2) of the doubled plasmid (fwd + rc)
       using Python str.find() — O(n), sub-millisecond for typical plasmids.
    2. Match on progressively shorter prefixes (50 → 30 → 15 AA) for
       robustness against occasional leading-residue differences.
    3. Falls back to Smith-Waterman only if no string match is found.

    Returns (gene_dna, gene_start, gene_length, is_rc_strand).
    gene_dna is always in 5'->3' coding orientation.
    gene_start / gene_length are coordinates on the original plus-strand plasmid.
    """
    plasmid = _clean(wt_plasmid_seq)
    protein = _clean(wt_protein_seq)
    L = len(plasmid)
    extended = plasmid + plasmid  # handle wrap-around

    for prefix_len in (min(50, len(protein)), min(30, len(protein)), min(15, len(protein))):
        prefix = protein[:prefix_len]
        for strand_seq, is_rc in [(extended, False), (_reverse_complement(extended), True)]:
            for frame in range(3):
                translated = _translate_full(strand_seq[frame:])
                pos = translated.find(prefix)
                if pos == -1:
                    continue

                # Expand to full protein length
                dna_start = frame + pos * 3
                dna_end   = frame + (pos + len(protein)) * 3

                if dna_end > len(strand_seq):
                    continue

                # gene_dna is sliced from strand_seq — already in coding orientation.
                gene_dna = strand_seq[dna_start:dna_end]

                # Map back to plus-strand coordinates on the original plasmid.
                # When is_rc=True, dna_start is a position in rc(extended).
                # The equivalent start on the plus-strand extended is (2L - dna_end).
                if is_rc:
                    gene_start = (2 * L - dna_end) % L
                else:
                    gene_start = dna_start % L
                gene_length = dna_end - dna_start

                gene_protein = _translate(gene_dna)
                print(f"  Gene found (prefix={prefix_len}AA, frame={frame}, rc={is_rc}): "
                      f"start={gene_start}, length={gene_length}bp ({len(gene_protein)} AA)")
                print(f"  First 50 AA: {gene_protein[:50]}")
                return gene_dna, gene_start, gene_length, is_rc

    # No exact prefix match -- fall back to Smith-Waterman
    print("  Fast search missed -- falling back to Smith-Waterman (this may take ~30s)...")
    return locate_gene_sw(wt_plasmid_seq, wt_protein_seq)


class SequenceAnalyzer:
    """
    Analyzes a batch of variant plasmids against a WT reference.
    Called once per experiment via analyze_variant_batch().

    Intentionally stateless — all data flows through method arguments and
    local variables.  The module-level singleton is therefore safe for
    concurrent calls from multiple background threads (one per experiment).
    """

    def analyze_variant_batch(
        self,
        variants_data: List[Dict],
        wt_protein_sequence: str,
        wt_plasmid_sequence: str,
    ) -> List[Dict]:
        """
        Locate WT gene once -- fast string search with SW fallback -- then extract
        the same DNA region from every variant and call identify_mutations().

        Each item in variants_data must contain:
            id, assembled_dna_sequence  (and any other keys to pass through)

        Returns the same list with 'protein_sequence' and 'mutations' added.
        """
        print(f"\nSequenceAnalyzer: processing {len(variants_data)} variants")
        print(f"  WT plasmid length : {len(wt_plasmid_sequence.strip())} bp")
        print(f"  WT protein length : {len(wt_protein_sequence.strip())} AA")

        # -- Step 1: locate gene in WT plasmid --
        # All gene-location data is kept in local variables so concurrent calls
        # for different experiments cannot overwrite each other's state.
        print("  Locating WT gene (fast search -> SW fallback)...")
        wt_gene_dna, gene_start, gene_length, is_gene_rc = locate_gene_fast(
            wt_plasmid_sequence, wt_protein_sequence
        )

        wt_protein = _translate(wt_gene_dna)
        plasmid_len = len(_clean(wt_plasmid_sequence))
        wraps = gene_start + gene_length > plasmid_len

        print(f"\n{'='*70}")
        print(f"LOCKED GENE COORDINATES")
        print(f"  start      : {gene_start}")
        print(f"  length     : {gene_length} bp  ({len(wt_protein)} AA)")
        print(f"  strand     : {'minus (RC)' if is_gene_rc else 'plus'}")
        print(f"  wraps orig : {wraps}")
        print(f"  plasmid len: {plasmid_len} bp")
        print(f"{'='*70}\n")

        # ── Step 2: process each variant ─────────────────────────────────────
        results = []
        for i, variant in enumerate(variants_data):
            if i % 50 == 0:
                print(f"  Processing variant {i+1}/{len(variants_data)}...")

            variant_copy = variant.copy()

            try:
                # Estimate circular rotation between WT and this variant assembly.
                # Variant plasmids may be sequenced/assembled starting at a different
                # position on the circle, shifting gene coordinates by a fixed offset.
                rotation = _estimate_rotation_offset(
                    wt_plasmid_sequence, variant['assembled_dna_sequence']
                )
                adj_start = (
                    (gene_start - rotation) % plasmid_len
                    if rotation is not None else gene_start
                )

                # extract_gene always returns plus-strand DNA for the region.
                # If the gene is on the minus strand, reverse-complement to get
                # the 5'->3' coding sequence before translating and comparing.
                var_gene_dna = extract_gene(
                    variant['assembled_dna_sequence'], adj_start, gene_length
                )
                if is_gene_rc:
                    var_gene_dna = _reverse_complement(var_gene_dna)

                var_protein = _translate(var_gene_dna)

                mutations = identify_mutations(wt_gene_dna, var_gene_dna)

                # Add aligned_position via Needleman-Wunsch global alignment.
                # When indels exist near or before a mutation site, simple 1-based
                # codon positions can be off by the indel count.  The aligned
                # position is used by the 3D fingerprint to correctly map mutations
                # onto AlphaFold structure residues.
                aligned_wt, aligned_var = _needleman_wunsch(wt_protein, var_protein)
                wt_pos_map = _build_wt_position_map(aligned_wt, aligned_var)
                for m in mutations:
                    m['aligned_position'] = wt_pos_map.get(m['position'], m['position'])

                # Debug first 3 variants
                if i < 3:
                    n_syn = sum(1 for m in mutations if m['mutation_type'] == 'synonymous')
                    n_ns = len(mutations) - n_syn
                    wt_len = len(wt_protein)
                    v_len = len(var_protein)
                    matches = sum(
                        1 for j in range(min(wt_len, v_len))
                        if wt_protein[j] == var_protein[j]
                    )
                    pct = matches / min(wt_len, v_len) * 100 if min(wt_len, v_len) else 0
                    print(f"  Variant {i+1}: {len(mutations)} mutations "
                          f"({n_ns} non-syn, {n_syn} syn), "
                          f"identity={pct:.1f}%  [{v_len} AA]")

                variant_copy['protein_sequence'] = var_protein
                variant_copy['mutations'] = mutations

            except Exception as e:
                print(f"  ERROR extracting variant {i+1}: {e}")
                variant_copy['protein_sequence'] = None
                variant_copy['mutations'] = []

            variant_copy['mutation_count'] = len(variant_copy['mutations'])
            results.append(variant_copy)

        # ── Step 3: summary ───────────────────────────────────────────────────
        total_syn = sum(
            sum(1 for m in v['mutations'] if m['mutation_type'] == 'synonymous')
            for v in results
        )
        total_ns = sum(
            sum(1 for m in v['mutations'] if m['mutation_type'] == 'non-synonymous')
            for v in results
        )
        avg_ns = total_ns / len(results) if results else 0
        print(f"\nAnalysis complete:")
        print(f"  Total non-synonymous : {total_ns}  (avg {avg_ns:.1f}/variant)")
        print(f"  Total synonymous     : {total_syn}")

        return results


# Module-level singleton imported by experiments.py.
# Safe for concurrent use — SequenceAnalyzer holds no mutable instance state.
sequence_analyzer = SequenceAnalyzer()
