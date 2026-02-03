"""
Staging service: fetch UniProt data and validate plasmid encodes WT protein.
"""
from __future__ import annotations

from dataclasses import asdict
from typing import Any, Optional

from .sequence_tools import parse_fasta_dna, parse_fasta_protein
from .plasmid_validation import find_wt_in_plasmid
from .uniprot_client import (
    UniProtError,
    fetch_uniprot_fasta,
    fetch_uniprot_features_json,
)


def stage_experiment_validate_plasmid(
    accession: str,
    plasmid_fasta_text: str,
    fetch_features: bool = True,
) -> dict[str, Any]:
    """
    Integration-ready “Part C staging” function.

    Web routes can call this directly later. CLI can call it too.
    Returns a JSON-friendly dict: {accession, wt_protein, features, validation, errors}.
    """
    result: dict[str, Any] = {
        "accession": accession.strip(),
        "wt_protein": None,
        "features": None,
        "validation": None,
        "error": None,
    }

    try:
        wt_fasta = fetch_uniprot_fasta(accession)
        wt_protein = parse_fasta_protein(wt_fasta)
        result["wt_protein"] = wt_protein

        features: Optional[list[dict[str, Any]]] = None
        if fetch_features:
            # Features may be empty; that’s okay. We store what UniProt provides.
            features = fetch_uniprot_features_json(accession)
        result["features"] = features

        plasmid_dna = parse_fasta_dna(plasmid_fasta_text)

        call = find_wt_in_plasmid(plasmid_dna, wt_protein)
        result["validation"] = asdict(call)

        return result

    except UniProtError as e:
        result["error"] = f"UniProt error: {e}"
        return result
    except Exception as e:
        result["error"] = f"Unexpected error: {e}"
        return result
