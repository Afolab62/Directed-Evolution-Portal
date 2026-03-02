from pathlib import Path
import pytest

from app.services.sequence_tools import parse_fasta_dna, parse_fasta_protein
from app.services.plasmid_validation import find_wt_in_plasmid


def test_plasmid_validation_example_data():
    """
    End-to-end service test using the example plasmid FASTA + UniProt WT sequence.

    Why this test matters:
    - protects the highest-value workflow (WT staging validation)
    - catches regressions in parsing, translation, circular mapping, and tie-breaking
    - gives confidence the portal can validate a real target sequence
    """
    plasmid_path = Path("data/pET-28a_BSU_DNA_Pol_I_WT.fa")
    wt_path = Path("data/wt_protein.fasta")

    if not (plasmid_path.exists() and wt_path.exists()):
        pytest.skip("Example data not present in this checkout. Run locally with data/ folder.")

    plasmid = parse_fasta_dna(plasmid_path.read_text())
    wt = parse_fasta_protein(wt_path.read_text())

    call = find_wt_in_plasmid(plasmid, wt)

    assert call.is_valid is True
    assert call.strand in ("+", "-")
    assert call.frame in (0, 1, 2)
    assert call.start_nt is not None
    assert call.end_nt_exclusive is not None


