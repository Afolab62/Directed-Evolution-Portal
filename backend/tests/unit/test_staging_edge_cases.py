from app.services.staging import stage_experiment_validate_plasmid
from app.services.uniprot_client import UniProtError


def test_staging_handles_invalid_fasta(monkeypatch):
    # Avoid real UniProt call
    import app.services.staging as staging

    def fake_fetch_uniprot_fasta(_accession: str, timeout_s: float = 10.0) -> str:
        return ">x\nMTEST"

    monkeypatch.setattr(staging, "fetch_uniprot_fasta", fake_fetch_uniprot_fasta)
    monkeypatch.setattr(staging, "fetch_uniprot_features_json", lambda _acc: [])

    bad_plasmid = ">p\nACGTXZZ"  # invalid DNA letters
    result = stage_experiment_validate_plasmid("O34996", bad_plasmid, fetch_features=False)

    assert result["validation"] is None
    assert result["error"] is not None


def test_staging_handles_uniprot_error(monkeypatch):
    import app.services.staging as staging

    def fake_fetch_uniprot_fasta(_accession: str, timeout_s: float = 10.0) -> str:
        raise UniProtError("Simulated UniProt failure")

    monkeypatch.setattr(staging, "fetch_uniprot_fasta", fake_fetch_uniprot_fasta)

    result = stage_experiment_validate_plasmid("BADACC", ">p\nAAA", fetch_features=False)

    assert result["validation"] is None
    assert result["error"] is not None
