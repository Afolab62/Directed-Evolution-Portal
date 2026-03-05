"""
mutation_analysis.py
--------------------
Part 4a — Sequence Analysis Layer
===================================
Pure sequence-analysis module: no Plotly, no Flask, no database.

Responsibilities
----------------
1. Locate the WT coding gene inside a (circular) plasmid.
2. Extract the equivalent gene region from any variant plasmid,
   accounting for circular rotation between assemblies.
3. Translate both gene sequences and compare them codon-by-codon
   to identify synonymous and non-synonymous mutations.
4. Perform a global protein alignment (Needleman-Wunsch, stdlib-only)
   so every mutation carries an alignment-aware position used by
   the visualisation layer.
5. Load variant tables from .tsv or .json files.
6. Export all mutations for a variant to a CSV file.

Public API (used by fingerprint_plot.py and the Flask app)
----------------------------------------------------------
    load_variants_table(path)          -> dict[int, dict]
    analyze_target_variant(...)        -> dict   (the "analysis" object)
    get_variant_lineage(...)           -> list[int]
    analyze_lineage_mutations(...)     -> list[dict]
    write_mutation_csv(path, analysis) -> None
    infer_uniprot_id(fasta_header)     -> str | None
    read_fasta(path)                   -> (header, sequence)
"""

from __future__ import annotations

import csv
import json
import logging
import re
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Genetic code
# ---------------------------------------------------------------------------

CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


# ---------------------------------------------------------------------------
# Sequence helpers
# ---------------------------------------------------------------------------

def clean_dna(seq: str) -> str:
    """Strip whitespace and keep only valid DNA bases."""
    return re.sub(r"[^ACGTN]", "", (seq or "").upper())


def clean_protein(seq: str) -> str:
    """Strip whitespace and keep only valid single-letter amino acid codes."""
    return re.sub(r"[^A-Z*]", "", (seq or "").upper())


def read_fasta(path: Path) -> tuple[str, str]:
    """Read a single-record FASTA file.  Returns (header_without_gt, sequence)."""
    lines = [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if not lines or not lines[0].startswith(">"):
        raise ValueError(f"Invalid FASTA format: {path}")
    return lines[0][1:], "".join(lines[1:])


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    trans = str.maketrans("ACGTN", "TGCAN")
    return clean_dna(seq).translate(trans)[::-1]


def translate_dna(dna_seq: str) -> str:
    """Translate a DNA sequence to protein (stops at first stop codon)."""
    dna = clean_dna(dna_seq)
    out: list[str] = []
    for i in range(0, len(dna) - 2, 3):
        aa = CODON_TABLE.get(dna[i: i + 3], "X")
        if aa == "*":
            break
        out.append(aa)
    return "".join(out)


# ---------------------------------------------------------------------------
# Gene location in a circular plasmid
# ---------------------------------------------------------------------------

def find_wt_gene_call(wt_plasmid_seq: str, wt_protein_seq: str) -> dict[str, Any]:
    """
    Locate the WT coding region inside a circular plasmid by performing an
    exact protein-sequence search across all 6 reading frames (both strands,
    3 frames each).  The plasmid is doubled before searching to handle genes
    that wrap around the origin.

    Returns a *gene_call* dict consumed by extract_gene_from_plasmid():
        strand          '+' or '-'
        frame           0 / 1 / 2
        start_nt        0-based nt start on the circular sequence
        end_nt_exclusive  0-based exclusive end
        wraps_origin    bool
        gene_nt_len     coding length in nt
    """
    plasmid = clean_dna(wt_plasmid_seq)
    wt_protein = clean_protein(wt_protein_seq).replace("*", "")

    if not plasmid or not wt_protein:
        raise ValueError("WT plasmid / protein sequence cannot be empty.")

    plasmid_len = len(plasmid)
    gene_nt_len = len(wt_protein) * 3
    dna2 = plasmid + plasmid  # doubled for circular search

    # Slide a codon-aligned window of gene_nt_len through each frame on each
    # strand.  Translating the full frame at once (the naive approach) fails
    # for complete plasmids because stop codons in other ORFs terminate the
    # translation before the gene of interest is reached.
    first_aa = wt_protein[0]
    candidates: list[tuple[int, int, str, int]] = []
    for strand in ("+", "-"):
        strand_seq = dna2 if strand == "+" else reverse_complement(dna2)
        max_start = len(strand_seq) - gene_nt_len
        for frame in (0, 1, 2):
            for nt_start in range(frame, max_start + 1, 3):
                # Fast first-residue filter before full translation
                if CODON_TABLE.get(strand_seq[nt_start:nt_start + 3], "X") != first_aa:
                    continue
                window = strand_seq[nt_start:nt_start + gene_nt_len]
                if translate_dna(window) != wt_protein:
                    continue
                if strand == "+":
                    start_nt2 = nt_start
                else:
                    start_nt2 = len(dna2) - (nt_start + gene_nt_len)
                rank = 0 if start_nt2 < plasmid_len else 1
                candidates.append((rank, start_nt2, strand, frame))

    if not candidates:
        raise ValueError(
            "Could not locate the WT protein sequence inside the WT plasmid. "
            "Verify that the correct plasmid and protein FASTA files were supplied."
        )

    candidates.sort(key=lambda x: (x[0], x[1], x[2], x[3]))
    _, start_nt2, strand, frame = candidates[0]
    start_nt = start_nt2 % plasmid_len
    wraps_origin = (start_nt2 < plasmid_len) and ((start_nt2 + gene_nt_len) > plasmid_len)

    if wraps_origin:
        end_nt_exclusive = (start_nt + gene_nt_len) - plasmid_len
    else:
        end_nt_exclusive = start_nt + gene_nt_len
        if end_nt_exclusive > plasmid_len:
            end_nt_exclusive %= plasmid_len

    logger.info(
        "WT gene located: strand=%s, frame=%d, start_nt=%d, "
        "gene_nt_len=%d, wraps_origin=%s",
        strand, frame, start_nt, gene_nt_len, wraps_origin,
    )

    return {
        "strand": strand,
        "frame": frame,
        "start_nt": start_nt,
        "end_nt_exclusive": end_nt_exclusive,
        "wraps_origin": wraps_origin,
        "gene_nt_len": gene_nt_len,
    }


def extract_gene_from_plasmid(plasmid_seq: str, gene_call: dict[str, Any]) -> str:
    """
    Slice the coding gene out of a plasmid using coordinates from gene_call.
    Handles origin-wrapping and reverse-strand genes.
    """
    plasmid = clean_dna(plasmid_seq)
    start_nt = int(gene_call["start_nt"])
    gene_nt_len = int(gene_call["gene_nt_len"])
    wraps_origin = bool(gene_call["wraps_origin"])
    strand = str(gene_call["strand"])

    if not wraps_origin:
        gene_seq = plasmid[start_nt: start_nt + gene_nt_len]
    else:
        overflow = gene_nt_len - (len(plasmid) - start_nt)
        gene_seq = plasmid[start_nt:] + plasmid[:overflow]

    if strand == "-":
        gene_seq = reverse_complement(gene_seq)

    return gene_seq


def extract_gene_with_start(
    plasmid_seq: str,
    start_nt: int,
    gene_nt_len: int,
    strand: str,
) -> str:
    """
    Variant of extract_gene_from_plasmid() accepting an explicit start offset.
    Used when a circular-rotation offset has been estimated between the WT
    and a variant plasmid assembly.
    """
    plasmid = clean_dna(plasmid_seq)
    if start_nt + gene_nt_len <= len(plasmid):
        gene_seq = plasmid[start_nt: start_nt + gene_nt_len]
    else:
        overflow = start_nt + gene_nt_len - len(plasmid)
        gene_seq = plasmid[start_nt:] + plasmid[:overflow]

    if strand == "-":
        gene_seq = reverse_complement(gene_seq)

    return gene_seq


# ---------------------------------------------------------------------------
# Circular-rotation estimation
# ---------------------------------------------------------------------------

def estimate_rotation_offset(wt_plasmid_seq: str, variant_plasmid_seq: str) -> int | None:
    """
    Estimate the circular rotation offset (in nt) between the WT and a
    variant plasmid assembly.  Multiple short anchors from the WT are located
    in the variant; the most-voted offset is returned.

    Returns None if no consistent offset can be determined.
    """
    wt = clean_dna(wt_plasmid_seq)
    var = clean_dna(variant_plasmid_seq)
    length = len(wt)

    anchor_positions = list(range(0, length, 350))
    anchor_sizes = (150, 120, 100, 80, 60, 40)

    for k in anchor_sizes:
        votes: dict[int, int] = {}
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


# ---------------------------------------------------------------------------
# Mutation detection
# ---------------------------------------------------------------------------

def analyze_mutations(wt_gene_seq: str, variant_gene_seq: str) -> list[dict[str, Any]]:
    """
    Codon-by-codon comparison of WT and variant gene sequences.

    Returns a list of mutation dicts, each containing:
        position        1-based codon position in the WT gene
        wt_codon        WT codon (3-letter DNA)
        mut_codon       variant codon
        wt_aa           WT amino acid (single letter)
        mut_aa          variant amino acid
        mutation_type   'synonymous' | 'non-synonymous'
        aa_change       human-readable label e.g. 'E51A'
    """
    wt = clean_dna(wt_gene_seq)
    var = clean_dna(variant_gene_seq)

    codon_len = min(len(wt), len(var))
    codon_len -= codon_len % 3

    out: list[dict[str, Any]] = []
    for i in range(0, codon_len, 3):
        wt_codon = wt[i: i + 3]
        mut_codon = var[i: i + 3]
        if wt_codon == mut_codon:
            continue

        pos = (i // 3) + 1
        wt_aa = CODON_TABLE.get(wt_codon, "X")
        mut_aa = CODON_TABLE.get(mut_codon, "X")
        mutation_type = "synonymous" if wt_aa == mut_aa else "non-synonymous"

        out.append({
            "position": pos,
            "wt_codon": wt_codon,
            "mut_codon": mut_codon,
            "wt_aa": wt_aa,
            "mut_aa": mut_aa,
            "mutation_type": mutation_type,
            "aa_change": f"{wt_aa}{pos}{mut_aa}",
        })

    return out


# ---------------------------------------------------------------------------
# Global protein alignment (Needleman-Wunsch, stdlib only)
# ---------------------------------------------------------------------------

def global_align(
    seq_a: str,
    seq_b: str,
    match: int = 2,
    mismatch: int = -1,
    gap: int = -2,
) -> tuple[str, str]:
    """
    Needleman-Wunsch global alignment of two protein sequences.
    Returns (aligned_seq_a, aligned_seq_b) with '-' gap characters.
    """
    n, m = len(seq_a), len(seq_b)

    score = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[0] * (m + 1) for _ in range(n + 1)]  # 0=diag, 1=up, 2=left

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
            s_diag = score[i - 1][j - 1] + (match if a == b else mismatch)
            s_up = score[i - 1][j] + gap
            s_left = score[i][j - 1] + gap

            best = s_diag
            t = 0
            if s_up > best:
                best = s_up
                t = 1
            if s_left > best:
                best = s_left
                t = 2
            score[i][j] = best
            trace[i][j] = t

    aligned_a: list[str] = []
    aligned_b: list[str] = []
    i, j = n, m

    while i > 0 or j > 0:
        t = trace[i][j] if i >= 0 and j >= 0 else 0
        if i > 0 and j > 0 and t == 0:
            aligned_a.append(seq_a[i - 1])
            aligned_b.append(seq_b[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or t == 1):
            aligned_a.append(seq_a[i - 1])
            aligned_b.append("-")
            i -= 1
        else:
            aligned_a.append("-")
            aligned_b.append(seq_b[j - 1])
            j -= 1

    aligned_a.reverse()
    aligned_b.reverse()
    return "".join(aligned_a), "".join(aligned_b)


def build_position_maps(
    aligned_wt: str,
    aligned_var: str,
) -> tuple[dict[int, int], dict[int, int]]:
    """
    Build position maps from an aligned pair.

    Returns:
        wt_map   {wt_residue_pos -> alignment_column}
        var_map  {var_residue_pos -> alignment_column}
    """
    wt_map: dict[int, int] = {}
    var_map: dict[int, int] = {}
    wt_pos = var_pos = 0
    for aln_pos, (wa, va) in enumerate(zip(aligned_wt, aligned_var), start=1):
        if wa != "-":
            wt_pos += 1
            wt_map[wt_pos] = aln_pos
        if va != "-":
            var_pos += 1
            var_map[var_pos] = aln_pos
    return wt_map, var_map


def alignment_stats(aligned_wt: str, aligned_var: str) -> dict[str, int]:
    """Count matches, mismatches, and gaps in an aligned pair."""
    matches = mismatches = gaps = 0
    for wa, va in zip(aligned_wt, aligned_var):
        if wa == "-" or va == "-":
            gaps += 1
        elif wa == va:
            matches += 1
        else:
            mismatches += 1
    return {"matches": matches, "mismatches": mismatches, "gaps": gaps}


# ---------------------------------------------------------------------------
# Variant table loading
# ---------------------------------------------------------------------------

def _parse_int(value: Any, field: str, default: int | None = None) -> int:
    if value is None or str(value).strip() == "":
        if default is not None:
            return default
        raise ValueError(f"Missing integer value for field '{field}'")
    s = str(value).strip()
    try:
        return int(float(s)) if "." in s else int(s)
    except Exception as exc:
        if default is not None:
            return default
        raise ValueError(f"Cannot parse integer for '{field}': {value!r}") from exc


def _pick_field(row: dict[str, Any], names: list[str], required: bool = True, default: Any = None) -> Any:
    for n in names:
        if n in row:
            return row[n]
    if required:
        raise KeyError(f"Required field not found. Tried: {', '.join(names)}")
    return default


def _normalize_variant_row(row: dict[str, Any]) -> dict[str, Any]:
    """Normalise a raw dict from the TSV/JSON into a consistent structure."""
    variant_id = _parse_int(
        _pick_field(row, ["Plasmid_Variant_Index", "Variant ID", "variant_id", "variantId"]),
        "variant_id",
    )
    parent_variant_id = _parse_int(
        _pick_field(row, ["Parent_Plasmid_Variant", "Parent", "parent_variant_id", "parent"],
                    required=False, default=-1),
        "parent_variant_id",
        default=-1,
    )
    generation = _parse_int(
        _pick_field(row, ["Directed_Evolution_Generation", "Generation", "generation"],
                    required=False, default=-1),
        "generation",
        default=-1,
    )
    dna_sequence = clean_dna(str(
        _pick_field(row, ["Assembled_DNA_Sequence", "DNA Sequence", "dna_sequence", "sequence"])
    ))

    if not dna_sequence:
        raise ValueError(f"Variant {variant_id} has an empty DNA sequence")

    protein_sequence = clean_protein(str(
        _pick_field(row, ["Protein_Sequence", "protein_sequence", "protein"], required=False, default="")
    )).replace("*", "")

    return {
        "variant_id": variant_id,
        "parent_variant_id": parent_variant_id,
        "generation": generation,
        "dna_sequence": dna_sequence,
        "protein_sequence": protein_sequence,
    }


def _load_rows_from_json(path: Path) -> list[dict[str, Any]]:
    obj = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(obj, list):
        return [r for r in obj if isinstance(r, dict)]
    if isinstance(obj, dict):
        for key in ("rows", "variants", "data"):
            if key in obj and isinstance(obj[key], list):
                return [r for r in obj[key] if isinstance(r, dict)]
        for value in obj.values():
            if isinstance(value, list) and value and isinstance(value[0], dict):
                return [r for r in value if isinstance(r, dict)]
    raise ValueError(f"Unsupported JSON structure in {path}")


def load_variants_table(path: Path) -> dict[int, dict[str, Any]]:
    """
    Load a variant table from a .tsv or .json file.

    Returns a dict keyed by integer variant ID.
    Each value is a normalised row dict with keys:
        variant_id, parent_variant_id, generation, dna_sequence
    """
    rows: list[dict[str, Any]] = []

    suffix = path.suffix.lower()
    if suffix in {".tsv", ".txt"}:
        with path.open("r", encoding="utf-8", newline="") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
    elif suffix == ".json":
        rows = _load_rows_from_json(path)
    else:
        raise ValueError(f"Unsupported file format '{suffix}'. Supply .tsv or .json")

    if not rows:
        raise ValueError(f"No rows loaded from {path}")

    variants: dict[int, dict[str, Any]] = {}
    for raw in rows:
        try:
            norm = _normalize_variant_row(raw)
            variants[norm["variant_id"]] = norm
        except (KeyError, ValueError) as exc:
            logger.warning("Skipping row — %s", exc)

    if not variants:
        raise ValueError(f"No valid variants found in {path}")

    logger.info("Loaded %d variants from %s", len(variants), path)
    return variants


# ---------------------------------------------------------------------------
# Main analysis entry point
# ---------------------------------------------------------------------------

def analyze_target_variant(
    variants: dict[int, dict[str, Any]],
    variant_id: int,
    wt_plasmid_seq: str,
    wt_protein_seq: str,
) -> dict[str, Any]:
    """
    Full analysis pipeline for a single variant.

    Steps
    -----
    1. Locate WT gene in WT plasmid (6-frame search).
    2. Estimate circular rotation between WT and variant plasmid.
    3. Extract variant gene using corrected coordinates.
    4. Translate both gene sequences.
    5. Detect mutations codon-by-codon.
    6. Global-align WT and variant proteins; add aligned_position to each mutation.
    7. Compute alignment statistics.

    Returns
    -------
    A dict containing all intermediate and final results needed by the
    visualisation layer (fingerprint_plot.py) and the Flask app.
    """
    if variant_id not in variants:
        sample_ids = sorted(variants.keys())[:12]
        raise ValueError(
            f"Variant ID {variant_id} not found in table. "
            f"Available IDs (first 12): {sample_ids}"
        )

    # Step 1 – locate WT gene
    gene_call = find_wt_gene_call(wt_plasmid_seq, wt_protein_seq)
    wt_gene_seq = extract_gene_from_plasmid(wt_plasmid_seq, gene_call)

    # Step 2 & 3 – extract variant gene (rotation-aware)
    target = variants[variant_id]
    variant_gene_seq, rotation_offset = _extract_variant_gene_sequence(
        target["dna_sequence"],
        wt_plasmid_seq,
        gene_call,
    )

    # Step 4 – translate
    wt_protein_clean = clean_protein(wt_protein_seq).replace("*", "")
    variant_protein_seq = translate_dna(variant_gene_seq)
    variant_protein_clean = clean_protein(variant_protein_seq).replace("*", "")

    # Step 5 – mutations
    mutations = analyze_mutations(wt_gene_seq, variant_gene_seq)

    # Step 6 – alignment & aligned positions
    aligned_wt, aligned_variant = global_align(wt_protein_clean, variant_protein_clean)
    wt_map, _ = build_position_maps(aligned_wt, aligned_variant)
    for m in mutations:
        m["aligned_position"] = wt_map.get(int(m["position"]), int(m["position"]))

    # Step 7 – alignment stats
    stats = alignment_stats(aligned_wt, aligned_variant)

    logger.info(
        "Variant %d: %d mutations (%d non-syn, %d syn)",
        variant_id,
        len(mutations),
        sum(1 for m in mutations if m["mutation_type"] != "synonymous"),
        sum(1 for m in mutations if m["mutation_type"] == "synonymous"),
    )

    return {
        # Identity
        "variant_id": variant_id,
        "parent_variant_id": int(target["parent_variant_id"]),
        "generation": int(target["generation"]),
        # Gene-level outputs
        "rotation_offset": rotation_offset,
        "gene_call": gene_call,
        "wt_gene_sequence": wt_gene_seq,
        "variant_gene_sequence": variant_gene_seq,
        # Protein-level outputs
        "wt_protein_sequence": wt_protein_clean,
        "variant_protein_sequence": variant_protein_clean,
        # Alignment
        "aligned_wt": aligned_wt,
        "aligned_variant": aligned_variant,
        "alignment_stats": stats,
        "alignment_length": len(aligned_wt),
        # Mutations (sorted by position)
        "mutations": sorted(mutations, key=lambda x: int(x["position"])),
        "num_mutations": len(mutations),
        "num_synonymous": sum(1 for m in mutations if m["mutation_type"] == "synonymous"),
        "num_nonsynonymous": sum(1 for m in mutations if m["mutation_type"] != "synonymous"),
    }


def _extract_variant_gene_sequence(
    variant_plasmid_seq: str,
    wt_plasmid_seq: str,
    gene_call: dict[str, Any],
) -> tuple[str, int | None]:
    """
    Extract the coding gene from a variant plasmid using the WT gene call and
    any inferred circular rotation offset.
    """
    rotation_offset = estimate_rotation_offset(wt_plasmid_seq, variant_plasmid_seq)

    if rotation_offset is None:
        return extract_gene_from_plasmid(variant_plasmid_seq, gene_call), None

    variant_gene_start = (
        int(gene_call["start_nt"]) - rotation_offset
    ) % len(clean_dna(wt_plasmid_seq))
    return (
        extract_gene_with_start(
            plasmid_seq=variant_plasmid_seq,
            start_nt=variant_gene_start,
            gene_nt_len=int(gene_call["gene_nt_len"]),
            strand=str(gene_call["strand"]),
        ),
        rotation_offset,
    )


def assign_generation_to_target_mutations(
    variants: dict[int, dict[str, Any]],
    variant_id: int,
    wt_plasmid_seq: str,
    wt_protein_seq: str,
    target_mutations: list[dict[str, Any]],
    default_generation: int | None = None,
) -> list[dict[str, Any]]:
    """
    Annotate the target variant's final codon-level mutations with the earliest
    lineage generation at which the final codon state appears.

    Unlike the protein-only lineage analysis, this preserves synonymous
    mutations so structure plots can render them.
    """
    fallback_generation = (
        int(default_generation)
        if default_generation is not None
        else int(variants.get(variant_id, {}).get("generation", -1))
    )

    if not target_mutations:
        return []

    lineage = get_variant_lineage(variants, variant_id)
    if not lineage:
        return [
            {
                **dict(mutation),
                "generation": fallback_generation,
                "variant_id": variant_id,
            }
            for mutation in target_mutations
        ]

    gene_call = find_wt_gene_call(wt_plasmid_seq, wt_protein_seq)
    prev_gene_seq = extract_gene_from_plasmid(wt_plasmid_seq, gene_call)
    target_by_position = {
        int(mutation["position"]): mutation
        for mutation in target_mutations
    }
    generation_by_position: dict[int, int] = {}

    for current_variant_id in lineage:
        current_variant = variants[current_variant_id]
        generation = int(current_variant["generation"])
        dna_sequence = str(current_variant.get("dna_sequence") or "")
        if not dna_sequence:
            logger.warning(
                "Variant %d has no dna_sequence; skipping codon lineage step.",
                current_variant_id,
            )
            continue

        try:
            current_gene_seq, _ = _extract_variant_gene_sequence(
                dna_sequence,
                wt_plasmid_seq,
                gene_call,
            )
        except Exception as exc:
            logger.warning(
                "Could not extract gene for variant %d; skipping codon lineage step (%s)",
                current_variant_id,
                exc,
            )
            continue

        if generation == 0:
            prev_gene_seq = current_gene_seq
            continue

        step_mutations = analyze_mutations(prev_gene_seq, current_gene_seq)
        for step_mutation in step_mutations:
            position = int(step_mutation["position"])
            final_mutation = target_by_position.get(position)
            if final_mutation is None:
                continue
            if step_mutation["mut_codon"] != final_mutation["mut_codon"]:
                continue
            generation_by_position.setdefault(position, generation)

        prev_gene_seq = current_gene_seq

    annotated_mutations: list[dict[str, Any]] = []
    for mutation in target_mutations:
        position = int(mutation["position"])
        annotated = dict(mutation)
        annotated["generation"] = generation_by_position.get(position, fallback_generation)
        annotated["variant_id"] = variant_id
        annotated_mutations.append(annotated)

    return sorted(
        annotated_mutations,
        key=lambda mutation: (int(mutation["generation"]), int(mutation["position"])),
    )


# ---------------------------------------------------------------------------
# CSV export
# ---------------------------------------------------------------------------

def write_mutation_csv(path: Path, analysis: dict[str, Any]) -> None:
    """
    Write all mutations for a variant to a CSV file.

    Columns: variant_id, generation, parent_variant_id,
             position, aligned_position, wt_aa, mut_aa,
             mutation_type, wt_codon, mut_codon, aa_change
    """
    fieldnames = [
        "variant_id", "generation", "parent_variant_id",
        "position", "aligned_position",
        "wt_aa", "mut_aa", "mutation_type",
        "wt_codon", "mut_codon", "aa_change",
    ]

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for m in analysis["mutations"]:
            writer.writerow({
                "variant_id": analysis["variant_id"],
                "generation": analysis["generation"],
                "parent_variant_id": analysis["parent_variant_id"],
                "position": m["position"],
                "aligned_position": m.get("aligned_position", m["position"]),
                "wt_aa": m["wt_aa"],
                "mut_aa": m["mut_aa"],
                "mutation_type": m["mutation_type"],
                "wt_codon": m["wt_codon"],
                "mut_codon": m["mut_codon"],
                "aa_change": m["aa_change"],
            })

    logger.info("Mutation CSV written: %s (%d rows)", path, len(analysis["mutations"]))


# ---------------------------------------------------------------------------
# Lineage analysis (used by the linear fingerprint plot)
# ---------------------------------------------------------------------------

def get_variant_lineage(
    variants: dict[int, dict[str, Any]],
    variant_id: int,
) -> list[int]:
    """
    Traverse parent_variant_id links from the target variant back to the root,
    then return the lineage in chronological order (root → target).

    The root is any variant whose parent_variant_id is -1 or not present in
    the table.  A cycle-guard prevents infinite loops in malformed tables.

    Returns a list of variant IDs, oldest ancestor first, e.g.
        [gen1_id, gen2_id, ..., variant_id]
    """
    lineage: list[int] = []
    current = variant_id
    seen: set[int] = set()

    while current != -1 and current in variants:
        if current in seen:
            logger.warning(
                "Cycle detected in variant lineage at ID %d; stopping traversal.",
                current,
            )
            break
        seen.add(current)
        lineage.append(current)
        current = int(variants[current]["parent_variant_id"])

    lineage.reverse()  # oldest ancestor first
    return lineage


def analyze_lineage_mutations(
    variants: dict[int, dict[str, Any]],
    variant_id: int,
    wt_plasmid_seq: str,
    wt_protein_seq: str,
) -> list[dict[str, Any]]:
    """
    Analyse mutations introduced at each step along the evolutionary lineage
    of the target variant.

    For each consecutive pair (parent → child) in the lineage, codon-level
    mutations are detected between the two gene sequences.  Each mutation dict
    is extended with:
        generation   int  — the child's generation number
        variant_id   int  — the child variant's ID

    The first variant in the lineage is always compared against the WT.

    Returns a flat list sorted by (generation, position).
    """
    lineage = get_variant_lineage(variants, variant_id)
    if not lineage:
        return []

    # Use the pre-computed protein sequences stored in the variants table.
    # Comparing protein sequences between consecutive lineage steps is reliable
    # and avoids the circular-plasmid rotation-offset issues that afflict raw
    # DNA comparison.  Synonymous mutations (same AA, different codon) are not
    # captured here, but they are invisible in an AA-level fingerprint anyway.
    wt_protein_clean = clean_protein(wt_protein_seq).replace("*", "")

    # Seed with WT; generation-0 variants are skipped below (they ARE the WT).
    prev_protein = wt_protein_clean

    all_mutations: list[dict[str, Any]] = []

    for vid in lineage:
        var = variants[vid]
        generation = int(var["generation"])

        # Generation 0 is the WT itself — use its protein as the new baseline
        # (in case it differs slightly from the FASTA file) then move on.
        if generation == 0:
            p = var.get("protein_sequence", "")
            if p:
                prev_protein = p
            continue

        curr_protein = var.get("protein_sequence", "")
        if not curr_protein:
            logger.warning(
                "Variant %d has no protein_sequence — skipping lineage step.", vid
            )
            continue

        # Align the two proteins before comparing so that insertions and
        # deletions don't cause a cascade of false positional mismatches.
        # A raw position-by-position comparison breaks whenever an indel
        # shifts all downstream residues, inflating mutation counts by
        # hundreds and suppressing text labels in the plot.
        aligned_prev, aligned_curr = global_align(prev_protein, curr_protein)

        prev_pos = 0  # 1-based position counter in prev_protein
        for ap, ac in zip(aligned_prev, aligned_curr):
            if ap != "-":
                prev_pos += 1

            if ap == "-":
                # Insertion in curr relative to prev
                all_mutations.append({
                    "position": prev_pos,
                    "aligned_position": prev_pos,
                    "wt_aa": "-",
                    "mut_aa": ac,
                    "mutation_type": "insertion",
                    "wt_codon": "n/a",
                    "mut_codon": "n/a",
                    "aa_change": f"->{ac}",
                    "generation": generation,
                    "variant_id": vid,
                })
            elif ac == "-":
                # Deletion in curr relative to prev
                all_mutations.append({
                    "position": prev_pos,
                    "aligned_position": prev_pos,
                    "wt_aa": ap,
                    "mut_aa": "-",
                    "mutation_type": "deletion",
                    "wt_codon": "n/a",
                    "mut_codon": "n/a",
                    "aa_change": f"{ap}{prev_pos}-",
                    "generation": generation,
                    "variant_id": vid,
                })
            elif ap != ac:
                # Substitution
                all_mutations.append({
                    "position": prev_pos,
                    "aligned_position": prev_pos,
                    "wt_aa": ap,
                    "mut_aa": ac,
                    "mutation_type": "non-synonymous",
                    "wt_codon": "n/a",
                    "mut_codon": "n/a",
                    "aa_change": f"{ap}{prev_pos}{ac}",
                    "generation": generation,
                    "variant_id": vid,
                })

        prev_protein = curr_protein

    logger.info(
        "Lineage for variant %d: %d step(s), %d total mutations",
        variant_id,
        len(lineage),
        len(all_mutations),
    )
    return sorted(all_mutations, key=lambda x: (x["generation"], x["position"]))


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------

def infer_uniprot_id(fasta_header: str) -> str | None:
    """
    Extract a UniProt accession from a FASTA header string.
    Handles both '|O34996|' pipe-delimited and bare accession formats.
    Returns None if no accession is found.
    """
    header = (fasta_header or "").strip()
    m = re.search(r"\|([A-Z0-9]{6,10})\|", header)
    if m:
        return m.group(1)
    token = header.split()[0].replace(">", "")
    if re.fullmatch(r"[A-Z0-9]{6,10}", token):
        return token
    return None