"""
Staging blueprint

Provides the HTTP interface for the experiment staging workflow:
  POST /staging/api/staging   →  fetch UniProt WT + validate plasmid
"""
from __future__ import annotations

from flask import Blueprint, request, jsonify

# Import so that tests can monkeypatch `staging_routes.stage_experiment_validate_plasmid`
from services.staging import stage_experiment_validate_plasmid  # noqa: F401

staging_bp = Blueprint("staging", __name__, url_prefix="/staging")


@staging_bp.post("/api/staging")
def staging_validate():
    """
    Validate that a plasmid encodes a UniProt WT protein.

    Request JSON body:
      accession      (str, required) — UniProt accession, e.g. "O34996"
      plasmid_fasta  (str, required) — FASTA text of the plasmid DNA
      fetch_features (bool, optional, default true)

    Response 200:
      {accession, wt_protein, features, validation, error}
    Response 400:
      {error: "..."}
    """
    # Hard cap: ~4 Mbp is far beyond any real expression vector (typical: 2–20 kbp).
    # Larger inputs risk timeouts in six-frame translation + fuzzy scan.
    _MAX_PLASMID_BYTES = 4_000_000

    body = request.get_json(silent=True) or {}

    accession: str = (body.get("accession") or "").strip()
    plasmid_fasta: str = (body.get("plasmid_fasta") or "").strip()
    fetch_features: bool = bool(body.get("fetch_features", True))

    if not accession or not plasmid_fasta:
        return jsonify({"error": "accession and plasmid_fasta are required"}), 400

    if len(plasmid_fasta) > _MAX_PLASMID_BYTES:
        return jsonify({
            "error": (
                f"Plasmid FASTA is too large ({len(plasmid_fasta):,} characters). "
                f"Maximum supported size is {_MAX_PLASMID_BYTES:,} characters "
                f"(\u2248{_MAX_PLASMID_BYTES // 1_000:,}\u202fkbp). "
                "Please ensure you have uploaded the correct file."
            )
        }), 413

    # Allow tests (and the route itself) to monkeypatch this module-level name
    import routes.staging as _self
    result = _self.stage_experiment_validate_plasmid(
        accession,
        plasmid_fasta,
        fetch_features=fetch_features,
    )

    return jsonify(result), 200
