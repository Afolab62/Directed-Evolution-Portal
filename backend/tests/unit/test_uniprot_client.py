from app.services.uniprot_client import UniProtNotFound, fetch_uniprot_fasta


class DummyHTTPError(Exception):
    def __init__(self, code: int):
        self.code = code


def test_uniprot_not_found(monkeypatch):
    # Monkeypatch urlopen inside the module to simulate a 404.
    import app.services.uniprot_client as uc

    def fake_urlopen(req, timeout=10.0):
        raise uc.HTTPError(req.full_url, 404, "Not Found", hdrs=None, fp=None)

    monkeypatch.setattr(uc, "urlopen", fake_urlopen)

    try:
        fetch_uniprot_fasta("BADACC")
        assert False, "Expected UniProtNotFound"
    except UniProtNotFound:
        assert True

