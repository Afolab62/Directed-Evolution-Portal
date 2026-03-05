from __future__ import annotations

import logging
from pathlib import Path

from flask import Blueprint, jsonify, render_template, request

from app.services.fingerprint_db import save_active_experiment
from app.services.staging import stage_experiment_validate_plasmid

staging_bp = Blueprint("staging", __name__, url_prefix="/staging")
logger = logging.getLogger(__name__)

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
_EXAMPLE_PLASMID_CANDIDATES = (
    _PROJECT_ROOT / "data" / "pET-28a_BSU_DNA_Pol_I_WT.fa",
    _PROJECT_ROOT / "pET-28a_BSU_DNA_Pol_I_WT.fa",
    _PROJECT_ROOT / "data" / "wt_dna.fasta",
)


def _get_plasmid_fasta_from_request() -> str:
    plasmid_fasta = request.form.get("plasmid_fasta", "").strip()
    uploaded = request.files.get("plasmid_file")
    if uploaded and uploaded.filename:
        plasmid_fasta = uploaded.read().decode("utf-8", errors="ignore")
    return plasmid_fasta


def _coerce_bool(value: object, *, default: bool = False) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _persist_result(result: dict) -> None:
    validation = result.get("validation") or {}
    if (
        result.get("error") is None
        and validation.get("is_valid")
        and result.get("accession")
        and result.get("wt_protein")
        and result.get("wt_plasmid_seq")
    ):
        save_active_experiment(
            uniprot_id=str(result["accession"]),
            wt_protein_seq=str(result["wt_protein"]),
            wt_plasmid_seq=str(result["wt_plasmid_seq"]),
        )
        logger.info("Saved staging result for accession=%s", result["accession"])


@staging_bp.get("/")
def staging_form():
    return render_template(
        "staging.html",
        result=None,
        accession="",
        plasmid_fasta="",
        fetch_features=True,
    )


@staging_bp.post("/")
def staging_submit():
    accession = request.form.get("accession", "").strip().upper()
    fetch_features = _coerce_bool(request.form.get("fetch_features"))
    plasmid_fasta = _get_plasmid_fasta_from_request()

    if not accession or not plasmid_fasta:
        return render_template(
            "staging.html",
            result={"error": "Provide a UniProt accession and plasmid FASTA (paste or upload)."},
            accession=accession,
            plasmid_fasta=plasmid_fasta,
            fetch_features=fetch_features,
        )

    result = stage_experiment_validate_plasmid(
        accession=accession,
        plasmid_fasta_text=plasmid_fasta,
        fetch_features=fetch_features,
    )
    _persist_result(result)

    logger.info(
        "staging_submit accession=%s valid=%s match_type=%s",
        accession,
        (result.get("validation") or {}).get("is_valid"),
        (result.get("validation") or {}).get("match_type"),
    )

    return render_template(
        "staging.html",
        result=result,
        accession=accession,
        plasmid_fasta=plasmid_fasta,
        fetch_features=fetch_features,
    )


@staging_bp.post("/api/staging")
def staging_api():
    data = request.get_json(silent=True) or {}
    accession = str(data.get("accession", "")).strip().upper()
    plasmid_fasta = str(data.get("plasmid_fasta", "")).strip()
    fetch_features = _coerce_bool(data.get("fetch_features"), default=True)

    if not accession or not plasmid_fasta:
        return jsonify({"error": "accession and plasmid_fasta are required"}), 400

    result = stage_experiment_validate_plasmid(
        accession=accession,
        plasmid_fasta_text=plasmid_fasta,
        fetch_features=fetch_features,
    )
    _persist_result(result)
    return jsonify(result), 200


@staging_bp.get("/example")
def staging_example():
    for path in _EXAMPLE_PLASMID_CANDIDATES:
        if path.exists():
            return jsonify({"fasta": path.read_text(encoding="utf-8")}), 200
    return jsonify({"error": "Example FASTA not found"}), 404