"""
app/routes/analysis.py
-----------------------
Flask routes for triggering and retrieving analysis results.

These routes are thin wrappers — they validate HTTP input,
call the service layer, and return JSON responses.
No SQL, no sequence logic here.
"""

from flask import Blueprint, request, jsonify, current_app
from app.services.analysis_service import run_variant_analysis, run_batch_analysis

analysis_bp = Blueprint("analysis", __name__, url_prefix="/analysis")


@analysis_bp.route("/variant/<string:variant_index>", methods=["POST"])
def analyse_single_variant(variant_index: str):
    """
    Trigger analysis for a single variant.

    POST /analysis/variant/<plasmid_variant_index>
    Body (JSON): { "experiment_id": 1 }

    Returns analysis results including mutation list and counts.
    """
    body = request.get_json(silent=True) or {}
    experiment_id = body.get("experiment_id")

    if experiment_id is None:
        return jsonify({"error": "experiment_id is required in the request body."}), 400

    db_path = current_app.config.get("DATABASE", "database.db")

    try:
        result = run_variant_analysis(
            experiment_id=int(experiment_id),
            plasmid_variant_index=variant_index,
            db_path=db_path,
        )
        # Convert MutationInfo dataclasses to dicts for JSON serialisation
        result["mutations"] = [vars(m) for m in result.get("mutations", [])]
        return jsonify(result), 200

    except ValueError as exc:
        return jsonify({"error": str(exc)}), 422

    except Exception as exc:
        current_app.logger.error(f"Unexpected error analysing variant {variant_index}: {exc}")
        return jsonify({"error": "Internal server error."}), 500


@analysis_bp.route("/batch", methods=["POST"])
def analyse_batch():
    """
    Trigger analysis for a batch of variants.

    POST /analysis/batch
    Body (JSON): {
        "experiment_id": 1,
        "variant_indices": ["VAR_001", "VAR_002", ...]
    }

    Returns a list of results. Failed variants include an 'error' key.
    """
    body = request.get_json(silent=True) or {}
    experiment_id = body.get("experiment_id")
    variant_indices = body.get("variant_indices", [])

    if experiment_id is None:
        return jsonify({"error": "experiment_id is required."}), 400

    if not variant_indices:
        return jsonify({"error": "variant_indices list is required and cannot be empty."}), 400

    db_path = current_app.config.get("DATABASE", "database.db")

    try:
        results = run_batch_analysis(
            experiment_id=int(experiment_id),
            variant_indices=variant_indices,
            db_path=db_path,
        )
        # Serialise MutationInfo dataclasses where present
        for r in results:
            if "mutations" in r:
                r["mutations"] = [vars(m) for m in r["mutations"]]
        return jsonify(results), 200

    except ValueError as exc:
        return jsonify({"error": str(exc)}), 422

    except Exception as exc:
        current_app.logger.error(f"Unexpected error in batch analysis: {exc}")
        return jsonify({"error": "Internal server error."}), 500
