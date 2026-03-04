"""
test_staging_api.py

Tests the staging HTTP endpoints using a minimal Flask test app that registers
only the staging blueprint — no database, no session setup required.
"""
from flask import Flask
import routes.staging as staging_routes
from routes.staging import staging_bp


def _make_app():
    """Minimal Flask app for testing — only the staging blueprint."""
    app = Flask(__name__)
    app.config["TESTING"] = True
    app.config["SECRET_KEY"] = "test-secret-key"
    app.register_blueprint(staging_bp)
    return app


def test_staging_api_success(monkeypatch):
    # Patch the staging call to avoid external API/network
    def fake_stage(accession: str, plasmid_fasta_text: str, fetch_features: bool = True):
        return {
            "accession": accession,
            "wt_protein": "FAKESEQ",
            "features": [],
            "validation": {"is_valid": True},
            "error": None,
        }

    monkeypatch.setattr(staging_routes, "stage_experiment_validate_plasmid", fake_stage)

    client = _make_app().test_client()

    resp = client.post(
        "/staging/api/staging",
        json={"accession": "O34996", "plasmid_fasta": ">p\nAAA", "fetch_features": False},
    )
    assert resp.status_code == 200
    data = resp.get_json()
    assert data["accession"] == "O34996"
    assert data["validation"]["is_valid"] is True
    assert data["error"] is None


def test_staging_api_missing_fields():
    client = _make_app().test_client()

    resp = client.post("/staging/api/staging", json={"accession": ""})
    assert resp.status_code == 400
