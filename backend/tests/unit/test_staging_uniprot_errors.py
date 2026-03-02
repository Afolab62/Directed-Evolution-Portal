from app.services.staging import stage_experiment_validate_plasmid


def test_staging_handles_uniprot_not_found(monkeypatch):
    """
    Staging should fail gracefully when UniProt accession is invalid / not found.

    We monkeypatch the UniProt fetch function so the test does not rely on network access.
    """
    import app.services.staging as staging

    def fake_fetch_uniprot_fasta(_accession: str, timeout_s: float = 10.0) -> str:
        raise Exception("Simulated UniProt failure (not found)")

    # Patch the symbol used inside staging.py (important: patch staging.fetch_uniprot_fasta)
    monkeypatch.setattr(staging, "fetch_uniprot_fasta", fake_fetch_uniprot_fasta)

    plasmid_fasta_text = ">p\n" + ("A" * 1000)

    result = stage_experiment_validate_plasmid("BADACC", plasmid_fasta_text, fetch_features=False)

    assert result["wt_protein"] is None
    assert result["validation"] is None
    assert result["error"] is not None
    assert "error" in result["error"].lower()

