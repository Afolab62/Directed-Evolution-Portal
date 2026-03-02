from app.services.plasmid_validation import find_wt_in_plasmid

CODON = {
    "A": "GCT",
    "C": "TGT",
    "D": "GAT",
    "E": "GAA",
    "F": "TTT",
    "G": "GGT",
    "H": "CAT",
    "I": "ATT",
    "K": "AAA",
    "L": "CTG",
    "M": "ATG",
    "N": "AAT",
    "P": "CCT",
    "Q": "CAA",
    "R": "CGT",
    "S": "TCT",
    "T": "ACT",
    "V": "GTG",
    "W": "TGG",
    "Y": "TAT",
}


def protein_to_dna(protein: str) -> str:
    return "".join(CODON[a] for a in protein)


def test_fuzzy_match_accepts_single_mismatch():
    # WT length >= 30 required by validator
    wt = "M" + "ACDEFGHIKLMNPQRSTVWY" + "ACDEFGHIKL"  # 31 aa
    wt_dna = protein_to_dna(wt)

    # Create a plasmid gene with ONE amino-acid substitution relative to WT.
    # Replace one codon to change an AA without changing length.
    # Here we mutate position 10 (0-based in protein) from original AA to 'A' codon GCT.
    mut_pos = 10
    mutated = wt_dna[: mut_pos * 3] + "GCT" + wt_dna[(mut_pos * 3) + 3 :]

    plasmid = "A" * 100 + mutated + "A" * 100

    call = find_wt_in_plasmid(plasmid, wt, fuzzy_fallback=True, min_identity=0.95)

    assert call.is_valid is True
    assert call.match_type in ("fuzzy", "exact")
    assert call.identity is not None
    assert call.identity >= 0.95


def test_exact_allows_x_wildcard_from_ambiguous_dna():
    wt = "M" + "ACDEFGHIKLMNPQRSTVWY" + "ACDEFGHIKL"  # 31 aa
    wt_dna = protein_to_dna(wt)

    # Introduce an ambiguous codon that translates to 'X' in the plasmid translation.
    # We replace one codon with 'NNN' so translation yields X at that position.
    ambig_pos = 5
    ambig = wt_dna[: ambig_pos * 3] + "NNN" + wt_dna[(ambig_pos * 3) + 3 :]

    plasmid = "A" * 100 + ambig + "A" * 100

    call = find_wt_in_plasmid(plasmid, wt, fuzzy_fallback=True)

    # Should still validate because translated 'X' is treated as wildcard.
    assert call.is_valid is True
    assert call.match_type == "exact"

