from __future__ import annotations

"""
UniProt client

What this module does:
- Fetch WT protein FASTA for a UniProt accession.
- (Optional) Fetch UniProt "features" annotations (JSON).
- Provide lightweight caching and clear, user-friendly error messages.
"""

from dataclasses import dataclass
from hashlib import sha1
from pathlib import Path
from typing import Any, Optional
import json
import time
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


# ----------------------------- Exceptions -----------------------------

class UniProtError(Exception):
    """Base error for UniProt client failures."""


class UniProtNotFound(UniProtError):
    """Accession invalid or not found (typically 400/404)."""


class UniProtNetworkError(UniProtError):
    """Network/HTTP errors that are not simple not-found."""


# ----------------------------- Data model -----------------------------

@dataclass(frozen=True)
class UniProtRecord:
    accession: str
    fasta_text: str
    features: Optional[list[dict[str, Any]]] = None


# ----------------------------- Constants -----------------------------

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"

# Disk cache lives under instance/ folder.
DEFAULT_CACHE_DIR = Path("instance") / "uniprot_cache"
DEFAULT_CACHE_TTL_S = 24 * 3600  # 24 hours


# ----------------------------- Cache helpers -----------------------------

def _cache_key(url: str, accept: str) -> str:
    h = sha1()
    h.update(url.encode("utf-8"))
    h.update(b"|")
    h.update(accept.encode("utf-8"))
    return h.hexdigest()


def _cache_paths(cache_dir: Path, key: str) -> tuple[Path, Path]:
    body_path = cache_dir / f"{key}.body"
    meta_path = cache_dir / f"{key}.meta.json"
    return body_path, meta_path


def _read_cache(cache_dir: Path, url: str, accept: str, ttl_s: float) -> Optional[str]:
    try:
        key = _cache_key(url, accept)
        body_path, meta_path = _cache_paths(cache_dir, key)
        if not body_path.exists() or not meta_path.exists():
            return None
        meta = json.loads(meta_path.read_text())
        ts = float(meta.get("ts", 0.0))
        if (time.time() - ts) > ttl_s:
            return None
        return body_path.read_text()
    except Exception:
        return None


def _write_cache(cache_dir: Path, url: str, accept: str, text: str) -> None:
    try:
        cache_dir.mkdir(parents=True, exist_ok=True)
        key = _cache_key(url, accept)
        body_path, meta_path = _cache_paths(cache_dir, key)
        body_path.write_text(text)
        meta_path.write_text(json.dumps({"ts": time.time(), "url": url, "accept": accept}))
    except Exception:
        return


# ----------------------------- HTTP core -----------------------------

def _http_get(
    url: str,
    *,
    accept: str,
    timeout_s: float,
    accession: Optional[str] = None,
    use_cache: bool = True,
    cache_dir: Path = DEFAULT_CACHE_DIR,
    cache_ttl_s: float = DEFAULT_CACHE_TTL_S,
) -> str:
    if use_cache:
        cached = _read_cache(cache_dir, url, accept, ttl_s=cache_ttl_s)
        if cached is not None:
            return cached

    req = Request(url, headers={"Accept": accept, "User-Agent": "Directed-Evolution-Portal/1.0"})

    try:
        with urlopen(req, timeout=timeout_s) as resp:
            raw = resp.read()
            text = raw.decode("utf-8", errors="replace")
    except HTTPError as e:
        if e.code in (400, 404):
            msg = "UniProt accession not found or invalid."
            if accession:
                msg = f"UniProt accession not found or invalid: {accession}"
            raise UniProtNotFound(msg) from e

        raise UniProtNetworkError(f"UniProt HTTP error {e.code}.") from e
    except (URLError, TimeoutError) as e:
        raise UniProtNetworkError("UniProt network error (timeout or connection failure).") from e

    if use_cache:
        _write_cache(cache_dir, url, accept, text)

    return text


# ----------------------------- Public API -----------------------------

def fetch_uniprot_fasta(accession: str, timeout_s: float = 10.0, use_cache: bool = True) -> str:
    """Fetch UniProt FASTA text for an accession."""
    acc = accession.strip()
    url = f"{UNIPROT_BASE}/{acc}.fasta"
    return _http_get(url, accept="text/plain", timeout_s=timeout_s, accession=acc, use_cache=use_cache)


def fetch_uniprot_features_json(accession: str, timeout_s: float = 10.0, use_cache: bool = True) -> list[dict[str, Any]]:
    """Fetch UniProt JSON and extract the 'features' list (if present)."""
    acc = accession.strip()
    url = f"{UNIPROT_BASE}/{acc}.json"
    text = _http_get(url, accept="application/json", timeout_s=timeout_s, accession=acc, use_cache=use_cache)
    data = json.loads(text)
    return list(data.get("features", []))


def fetch_uniprot_protein_metadata(accession: str, timeout_s: float = 10.0, use_cache: bool = True) -> dict[str, Any]:
    """Fetch UniProt JSON and return full protein metadata."""
    acc = accession.strip()
    url = f"{UNIPROT_BASE}/{acc}.json"
    text = _http_get(url, accept="application/json", timeout_s=timeout_s, accession=acc, use_cache=use_cache)
    return json.loads(text)


def fetch_uniprot_record(accession: str, *, fetch_features: bool = True, timeout_s: float = 10.0,
                         use_cache: bool = True) -> UniProtRecord:
    """Convenience wrapper to fetch FASTA and (optionally) features together."""
    fasta = fetch_uniprot_fasta(accession, timeout_s=timeout_s, use_cache=use_cache)
    features = fetch_uniprot_features_json(accession, timeout_s=timeout_s, use_cache=use_cache) if fetch_features else None
    return UniProtRecord(accession=accession.strip(), fasta_text=fasta, features=features)


# ----------------------------- Extended / detailed fetch ---------------------

def _extract_go_terms(data: dict[str, Any]) -> list[dict[str, str]]:
    """Extract GO term annotations from a UniProt JSON record."""
    go_terms: list[dict[str, str]] = []
    aspect_map = {"F": "molecular_function", "C": "cellular_component", "P": "biological_process"}
    for ref in data.get("uniProtKBCrossReferences", []):
        if ref.get("database") != "GO":
            continue
        go_id = ref.get("id", "")
        term_name = ""
        aspect = ""
        for prop in ref.get("properties", []):
            if prop.get("key") == "GoTerm":
                raw = prop.get("value", "")
                if len(raw) >= 2 and raw[1] == ":":
                    aspect = aspect_map.get(raw[0], raw[0])
                    term_name = raw[2:]
                else:
                    term_name = raw
        if go_id:
            go_terms.append({"id": go_id, "term": term_name, "aspect": aspect})
    return go_terms


_DETAIL_DATABASES = frozenset((
    "PDB", "InterPro", "Pfam", "AlphaFoldDB", "KEGG",
    "RefSeq", "Ensembl", "PANTHER", "Reactome",
))


def fetch_uniprot_detailed(
    accession: str,
    timeout_s: float = 15.0,
    use_cache: bool = True,
) -> dict[str, Any]:
    """
    Return enriched protein data including cross-references, GO terms, PDB
    structure IDs, AlphaFold ID, InterPro/Pfam domains, KEGG IDs, keywords,
    taxonomy lineage, and gene names.  Built on top of the standard JSON
    fetch — no extra HTTP request is needed.
    """
    data = fetch_uniprot_protein_metadata(accession, timeout_s=timeout_s, use_cache=use_cache)

    # --- Cross-references grouped by database ---
    crossrefs_raw: dict[str, list[str]] = {}
    for ref in data.get("uniProtKBCrossReferences", []):
        db = ref.get("database", "")
        ref_id = ref.get("id", "")
        if db and ref_id:
            crossrefs_raw.setdefault(db, []).append(ref_id)

    crossrefs = {db: ids for db, ids in crossrefs_raw.items() if db in _DETAIL_DATABASES}

    # --- GO terms ---
    go_terms = _extract_go_terms(data)

    # --- Keywords ---
    keywords = [kw.get("name", "") for kw in data.get("keywords", []) if kw.get("name")]

    # --- Gene names (primary + synonyms) ---
    gene_names: list[str] = []
    for gene in data.get("genes", []):
        if "geneName" in gene:
            gene_names.append(gene["geneName"].get("value", ""))
        for syn in gene.get("synonyms", []):
            val = syn.get("value", "")
            if val:
                gene_names.append(val)

    # --- Taxonomy lineage ---
    # UniProt REST API returns lineage as a list of strings
    lineage = [
        (t if isinstance(t, str) else t.get("scientificName", ""))
        for t in data.get("organism", {}).get("lineage", [])
        if t
    ]

    return {
        "pdb_ids":           crossrefs.get("PDB", []),
        "alphafold_id":      (crossrefs.get("AlphaFoldDB") or [None])[0],
        "interpro_ids":      crossrefs.get("InterPro", []),
        "pfam_ids":          crossrefs.get("Pfam", []),
        "kegg_ids":          crossrefs.get("KEGG", []),
        "go_terms":          go_terms,
        "keywords":          keywords,
        "gene_names":        gene_names,
        "taxonomy_lineage":  lineage,
        "crossrefs":         crossrefs,
    }
