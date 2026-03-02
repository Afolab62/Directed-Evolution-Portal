from app.services.plasmid_validation import find_wt_in_plasmid

CODON = {
    "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT",
    "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAA", "L": "CTG",
    "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT",
    "S": "TCT", "T": "ACT", "V": "GTG", "W": "TGG", "Y": "TAT",
}

def protein_to_dna(protein: str) -> str:
    return "".join(CODON[a] for a in protein)

def test_wraparound_detection():
    wt = "M" + "ACDEFGHIKLMNPQRSTVWY"  # 21 aa
    gene = protein_to_dna(wt)

    # Split gene so it wraps around the end/start of the plasmid
    split = 10 * 3
    plasmid = gene[split:] + ("A" * 50) + gene[:split]

    # Lower min_wt_len for short synthetic WT sequence used in this test.
    call = find_wt_in_plasmid(plasmid, wt, min_wt_len=10)

    assert call.is_valid is True
    assert call.wraps_origin is True
