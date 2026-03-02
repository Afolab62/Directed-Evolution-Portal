from app import create_app


def test_staging_api_success(monkeypatch):
    # Patch the staging call to avoid external API/network
    import app.routes.staging as staging_routes

    def fake_stage(accession: str, plasmid_fasta_text: str, fetch_features: bool = True):
        return {
            "accession": accession,
            "wt_protein": "FAKESEQ",
            "features": [],
            "validation": {"is_valid": True},
            "error": None,
        }

    monkeypatch.setattr(staging_routes, "stage_experiment_validate_plasmid", fake_stage)

    app = create_app()
    app.config["TESTING"] = True
    client = app.test_client()

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
    app = create_app()
    app.config["TESTING"] = True
    client = app.test_client()

    resp = client.post("/staging/api/staging", json={"accession": ""})
    assert resp.status_code == 400
