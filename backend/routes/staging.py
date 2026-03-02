"""
Staging routes: UniProt accession + plasmid FASTA upload + validation.
"""
from pathlib import Path
import logging

from flask import Blueprint, jsonify, request

from services.staging import stage_experiment_validate_plasmid

staging_bp = Blueprint("staging", __name__, url_prefix="/staging")
logger = logging.getLogger(__name__)


def _get_plasmid_fasta_from_request() -> str:
    """Extract plasmid FASTA text from multipart file upload or form field."""
    plasmid_fasta = request.form.get("plasmid_fasta", "").strip()
    uploaded = request.files.get("plasmid_file")
    if uploaded and uploaded.filename:
        plasmid_fasta = uploaded.read().decode("utf-8")
    return plasmid_fasta


@staging_bp.post("/")
def staging_submit():
    """
    Form-data or multipart submission.
    Accepts: accession (str), plasmid_fasta (str) or plasmid_file (file), fetch_features (bool).
    Returns JSON validation result.
    """
    accession = request.form.get("accession", "").strip()
    fetch_features = bool(request.form.get("fetch_features", True))
    plasmid_fasta = _get_plasmid_fasta_from_request()

    if not accession or not plasmid_fasta:
        return jsonify({"error": "accession and plasmid_fasta are required"}), 400

    result = stage_experiment_validate_plasmid(
        accession=accession,
        plasmid_fasta_text=plasmid_fasta,
        fetch_features=fetch_features,
    )

    logger.info(
        "staging_submit accession=%s valid=%s match_type=%s",
        accession,
        (result.get("validation") or {}).get("is_valid"),
        (result.get("validation") or {}).get("match_type"),
    )

    return jsonify(result), 200


@staging_bp.post("/api/staging")
def staging_api():
    """
    JSON API endpoint for staging.
    Expects: {accession, plasmid_fasta, fetch_features?}
    Returns: {accession, wt_protein, wt_plasmid_seq, features, validation, error}
    """
    data = request.get_json(silent=True) or {}
    accession = str(data.get("accession", "")).strip()
    plasmid_fasta = str(data.get("plasmid_fasta", "")).strip()
    fetch_features = bool(data.get("fetch_features", True))

    if not accession or not plasmid_fasta:
        return jsonify({"error": "accession and plasmid_fasta are required"}), 400

    result = stage_experiment_validate_plasmid(
        accession=accession,
        plasmid_fasta_text=plasmid_fasta,
        fetch_features=fetch_features,
    )

    logger.info(
        "staging_api accession=%s valid=%s match_type=%s",
        accession,
        (result.get("validation") or {}).get("is_valid"),
        (result.get("validation") or {}).get("match_type"),
    )

    return jsonify(result), 200


@staging_bp.get("/example")
def staging_example():
    """Return example plasmid FASTA text for quick demos."""
    example_path = Path("data") / "pET-28a_BSU_DNA_Pol_I_WT.fa"
    if not example_path.exists():
        return jsonify({"error": "Example FASTA not found"}), 404
    return jsonify({"fasta": example_path.read_text(encoding="utf-8")}), 200
