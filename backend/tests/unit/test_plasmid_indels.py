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

def test_alignment_fallback_handles_insertion():
    wt = "M" + "ACDEFGHIKLMNPQRSTVWY" + "ACDEFGHIKL"  # 31 aa
    gene = protein_to_dna(wt)

    # Insert 1 codon (3 nt) inside the gene -> causes an AA insertion relative to WT
    ins_pos = 12
    gene_ins = gene[: ins_pos * 3] + "GCT" + gene[ins_pos * 3 :]

    plasmid = "A" * 100 + gene_ins + "A" * 100

    call = find_wt_in_plasmid(plasmid, wt)

    assert call.is_valid is True
    assert call.match_type in ("align", "exact", "fuzzy")
    # Exact/fuzzy should normally fail because lengths differ; align should succeed.
    assert call.identity is None or call.identity >= 0.90

def test_alignment_fallback_handles_deletion():
    wt = "M" + "ACDEFGHIKLMNPQRSTVWY" + "ACDEFGHIKL"  # 31 aa
    gene = protein_to_dna(wt)

    # Delete 1 codon (3 nt) inside the gene -> causes an AA deletion relative to WT
    del_pos = 8
    gene_del = gene[: del_pos * 3] + gene[(del_pos * 3) + 3 :]

    plasmid = "A" * 100 + gene_del + "A" * 100

    call = find_wt_in_plasmid(plasmid, wt)

    assert call.is_valid is True
    assert call.match_type in ("align", "exact", "fuzzy")
    assert call.identity is None or call.identity >= 0.90

