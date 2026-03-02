from app.services.plasmid_validation import find_wt_in_plasmid


# Simple codon choices to construct a predictable coding sequence (no stops).
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


def test_wraparound_match_sets_wraps_origin_true():
    # Ensure protein length >= 30 because the validator guards against short sequences.
    wt_protein = "M" + "ACDEFGHIKLMNPQRSTVWY" + "ACDEFGHIKL"  # 1 + 20 + 10 = 31 aa
    gene_dna = protein_to_dna(wt_protein)

    # Build a circular plasmid where the gene crosses the end/start boundary:
    # put the last part of gene at the end, first part at the beginning.
    L = 600  # plasmid length
    filler = "A" * (L - len(gene_dna))  # safe filler in frame +0
    split = 45  # nt split point inside gene_dna (must be multiple of 3 for frame stability)
    split -= (split % 3)

    end_chunk = gene_dna[split:]    # goes at end
    start_chunk = gene_dna[:split]  # goes at start

    plasmid = end_chunk + filler + start_chunk

    call = find_wt_in_plasmid(plasmid, wt_protein)

    assert call.is_valid is True
    assert call.wraps_origin is True
    assert call.strand == "+"
    assert call.frame == 0


def test_no_match_returns_is_valid_false():
    plasmid = "A" * 2000
    wt_protein = "M" + "ACDEFGHIKLMNPQRSTVWY" + "ACDEFGHIKL"  # 31 aa
    call = find_wt_in_plasmid(plasmid, wt_protein)

    assert call.is_valid is False
    assert call.match_type == "none"

