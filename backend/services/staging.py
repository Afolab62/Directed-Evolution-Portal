"""
Staging service: fetch UniProt data and validate plasmid encodes WT protein.
"""
from __future__ import annotations

from dataclasses import asdict
from typing import Any, Optional

from .sequence_tools import parse_fasta_dna, parse_fasta_protein
from .errors import FastaParseError, InvalidSequenceError
from .plasmid_validation import find_wt_in_plasmid
from .uniprot_client import (
    UniProtError,
    fetch_uniprot_fasta,
    fetch_uniprot_protein_metadata,
)


def stage_experiment_validate_plasmid(
    accession: str,
    plasmid_fasta_text: str,
    fetch_features: bool = True,
) -> dict[str, Any]:
    """
    Integration-ready staging function.

    Returns a JSON-friendly dict: {accession, wt_protein, features, validation, error}.
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

        # Fetch full UniProt record with protein metadata
        features_data: Optional[dict[str, Any]] = None
        if fetch_features:
            protein_metadata = fetch_uniprot_protein_metadata(accession)
            # Extract relevant data from the full metadata
            features_data = {
                "name": protein_metadata.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value"),
                "organism": protein_metadata.get("organism", {}).get("scientificName"),
                "features": protein_metadata.get("features", [])
            }
        result["features"] = features_data

        plasmid_dna = parse_fasta_dna(plasmid_fasta_text)

        call = find_wt_in_plasmid(plasmid_dna, wt_protein)
        result["validation"] = asdict(call)

        return result

    except UniProtError as e:
        result["error"] = f"UniProt error: {e}"
        return result
    except FastaParseError as e:
        result["error"] = (
            f"Invalid plasmid file: {e}. "
            "Please upload a FASTA file containing a single DNA sequence "
            "(must begin with a '>' header line)."
        )
        return result
    except InvalidSequenceError as e:
        result["error"] = (
            f"Plasmid sequence contains invalid characters: {e}. "
            "DNA sequences may only contain A, C, G, T and N (ambiguous base)."
        )
        return result
    except Exception as e:
        result["error"] = f"Unexpected error: {e}"
        return result
