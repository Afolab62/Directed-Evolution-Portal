"""
services/fingerprint_plot.py
----------------------------
3D mutation fingerprint figure builder.

Downloads and caches an AlphaFold/AlphaFill structure for the experiment
protein, parses Cα coordinates with pLDDT confidence scores, and builds a
Plotly 3D figure that maps variant mutations onto the protein backbone.

Public API
----------
    resolve_structure(uniprot_id, no_3d=False)
        -> (backbone, pdb_coords, info)
           backbone   : list[{r, x, y, z, plddt}]   — sorted Cα trace
           pdb_coords : dict[int, {x, y, z, plddt}]  — per-residue lookup
           info       : {status, source}

    build_3d_fingerprint(lineage_mutations, backbone, pdb_coords, variant_id,
                         uniprot_id, structure_source, feature_annotations)
        -> plotly.Figure
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import requests
import plotly.graph_objects as go

logger = logging.getLogger(__name__)

# Persistent PDB/CIF cache — survives backend restarts
_CACHE_DIR = (
    Path(__file__).resolve().parent.parent.parent / "instance" / "pdb_cache"
)



# ---------------------------------------------------------------------------
# Structure download helpers
# ---------------------------------------------------------------------------

def _get(url: str, timeout: float = 30.0) -> str:
    r = requests.get(url, timeout=timeout)
    r.raise_for_status()
    return r.text


def _fetch_alphafold_pdb(uniprot_id: str) -> tuple[Path, str]:
    """Download from AlphaFold EBI — PDB with pLDDT in B-factor column."""
    _CACHE_DIR.mkdir(parents=True, exist_ok=True)
    for version in ("v4", "v3", "v2", "v1"):
        fname = f"AF-{uniprot_id}-F1-model_{version}.pdb"
        cached = _CACHE_DIR / fname
        if cached.exists():
            return cached, fname
        url = f"https://alphafold.ebi.ac.uk/files/{fname}"
        try:
            text = _get(url)
        except Exception:
            continue
        if "ATOM" not in text:
            continue
        cached.write_text(text, encoding="utf-8")
        return cached, fname
    raise FileNotFoundError(f"No AlphaFold EBI model for '{uniprot_id}'")


def _fetch_alphafill(uniprot_id: str) -> tuple[Path, str]:
    """Download from AlphaFill — mmCIF with B_iso_or_equiv = pLDDT."""
    _CACHE_DIR.mkdir(parents=True, exist_ok=True)
    for version in ("v4", "v3", "v2", "v1"):
        model_id = f"AF-{uniprot_id}-F1-model_{version}"
        fname = f"{model_id}.cif"
        cached = _CACHE_DIR / fname
        if cached.exists():
            return cached, f"AlphaFill ({model_id})"
        url = f"https://alphafill.eu/v1/aff/{model_id}"
        try:
            text = _get(url)
        except Exception:
            continue
        if "_atom_site." not in text:
            continue
        cached.write_text(text, encoding="utf-8")
        return cached, f"AlphaFill ({model_id})"
    raise FileNotFoundError(f"No AlphaFill model for '{uniprot_id}'")


def _fetch_uniprot_pdb(uniprot_id: str) -> tuple[Path, str]:
    """Fetch the first PDB cross-referenced in UniProt."""
    _CACHE_DIR.mkdir(parents=True, exist_ok=True)
    data = requests.get(
        f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json", timeout=20
    ).json()
    pdb_ids = [
        ref["id"]
        for ref in data.get("uniProtKBCrossReferences", [])
        if ref.get("database") == "PDB" and ref.get("id")
    ]
    for pdb_id in pdb_ids:
        cached = _CACHE_DIR / f"{pdb_id}.pdb"
        if cached.exists():
            return cached, f"PDB ({pdb_id})"
        try:
            text = _get(f"https://files.rcsb.org/download/{pdb_id}.pdb")
        except Exception:
            continue
        if "ATOM" not in text:
            continue
        cached.write_text(text, encoding="utf-8")
        return cached, f"PDB ({pdb_id})"
    raise FileNotFoundError(f"No PDB entry for '{uniprot_id}'")


# ---------------------------------------------------------------------------
# Structure parsers — return (backbone, pdb_coords) with pLDDT
# ---------------------------------------------------------------------------

def _parse_pdb(path: Path) -> tuple[list[dict], dict[int, dict]]:
    """Parse PDB Cα atoms; pLDDT from B-factor column 60-66."""
    coords: dict[int, dict] = {}
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[12:16].strip() != "CA":
            continue
        try:
            res = int(line[22:26])
            coords[res] = {
                "x": round(float(line[30:38]), 3),
                "y": round(float(line[38:46]), 3),
                "z": round(float(line[46:54]), 3),
                "plddt": round(float(line[60:66]), 1),
            }
        except (ValueError, IndexError):
            continue
    if not coords:
        raise ValueError(f"No Cα atoms in {path}")
    backbone = [{"r": r, **c} for r, c in sorted(coords.items())]
    return backbone, coords


def _parse_mmcif(path: Path) -> tuple[list[dict], dict[int, dict]]:
    """Parse mmCIF Cα atoms; pLDDT from B_iso_or_equiv column."""
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    atom_cols: list[str] = []
    data_start: int | None = None

    for i, line in enumerate(lines):
        if line.strip() != "loop_":
            continue
        cols: list[str] = []
        j = i + 1
        while j < len(lines) and lines[j].lstrip().startswith("_"):
            cols.append(lines[j].strip())
            j += 1
        if cols and cols[0].startswith("_atom_site."):
            atom_cols = cols
            data_start = j
            break

    if not atom_cols or data_start is None:
        raise ValueError(f"No atom_site loop in {path}")

    col_idx = {name: i for i, name in enumerate(atom_cols)}

    def _pick(*names: str) -> int:
        for n in names:
            if n in col_idx:
                return col_idx[n]
        raise KeyError(names[0])

    atom_i  = _pick("_atom_site.auth_atom_id",  "_atom_site.label_atom_id")
    chain_i = _pick("_atom_site.auth_asym_id",  "_atom_site.label_asym_id")
    seq_i   = _pick("_atom_site.auth_seq_id",   "_atom_site.label_seq_id")
    x_i     = _pick("_atom_site.Cartn_x")
    y_i     = _pick("_atom_site.Cartn_y")
    z_i     = _pick("_atom_site.Cartn_z")
    b_i     = col_idx.get("_atom_site.B_iso_or_equiv")

    chains: dict[str, dict[int, dict]] = {}
    for line in lines[data_start:]:
        s = line.strip()
        if not s or s == "#" or s == "loop_" or s.startswith("_"):
            if chains:
                break
            continue
        parts = s.split()
        if len(parts) < len(atom_cols):
            continue
        try:
            if parts[atom_i] != "CA":
                continue
            chain_id = parts[chain_i] if parts[chain_i] not in {".", "?"} else "_"
            seq_s = parts[seq_i]
            if seq_s in {".", "?"}:
                continue
            pos    = int(seq_s)
            x      = float(parts[x_i])
            y      = float(parts[y_i])
            z      = float(parts[z_i])
            plddt  = round(float(parts[b_i]), 1) if b_i is not None else 70.0
        except (ValueError, IndexError):
            continue
        chains.setdefault(chain_id, {})[pos] = {
            "x": round(x, 3), "y": round(y, 3), "z": round(z, 3), "plddt": plddt
        }

    if not chains:
        raise ValueError(f"No Cα atoms in {path}")

    _, best = max(chains.items(), key=lambda kv: len(kv[1]))
    backbone = [{"r": r, **c} for r, c in sorted(best.items())]
    return backbone, best


def _parse_structure(path: Path) -> tuple[list[dict], dict[int, dict]]:
    """Auto-detect PDB vs mmCIF and parse."""
    text = path.read_text(encoding="utf-8", errors="ignore")[:4096]
    if "_atom_site." in text:
        return _parse_mmcif(path)
    return _parse_pdb(path)


# ---------------------------------------------------------------------------
# Public: resolve_structure
# ---------------------------------------------------------------------------

def resolve_structure(
    uniprot_id: str | None,
    no_3d: bool = False,
) -> tuple[list[dict], dict[int, dict], dict[str, str]]:
    """
    Download (or load from cache) a protein structure and parse Cα atoms.

    Priority: AlphaFold EBI PDB → AlphaFill mmCIF → UniProt cross-ref PDB

    Returns
    -------
    backbone    list[{r, x, y, z, plddt}]   — sorted Cα trace
    pdb_coords  dict[int, {x, y, z, plddt}] — per-residue lookup
    info        {status, source}
    """
    info: dict[str, str] = {}

    if no_3d:
        info["status"] = "3D disabled"
        return [], {}, info

    if not uniprot_id:
        info["status"] = "No UniProt ID — 3D unavailable"
        return [], {}, info

    for fetch_fn, label in (
        (_fetch_alphafold_pdb, "AlphaFold EBI"),
        (_fetch_alphafill,     "AlphaFill"),
        (_fetch_uniprot_pdb,   "UniProt/PDB"),
    ):
        try:
            path, source = fetch_fn(uniprot_id)
            backbone, pdb_coords = _parse_structure(path)
            info.update(status="3D mapping enabled", source=source)
            return backbone, pdb_coords, info
        except Exception as exc:
            logger.warning("%s failed for %s: %s", label, uniprot_id, exc)

    info["status"] = "3D unavailable (download failed)"
    return [], {}, info


# ---------------------------------------------------------------------------
# Functional region overlay helpers
# ---------------------------------------------------------------------------

_FUNC_COLORS: dict[str, str] = {
    "Exonuclease":     "#14b8a6",
    "Polymerase core": "#f59e0b",
    "Catalytic site":  "#f43f5e",
}


def _extract_functional_annotations(
    annotations: list[dict] | None,
) -> tuple[list[dict], list[dict]]:
    if not annotations:
        return [], []

    regions: list[dict] = []
    sites:   list[dict] = []
    seen_r: set[tuple] = set()
    seen_s: set[tuple] = set()

    for feat in annotations:
        ftype = str(feat.get("type") or "").strip().lower()
        desc  = str(feat.get("description") or "").strip()
        text  = f"{ftype} {desc}".lower()
        start = feat.get("start")
        end   = feat.get("end")
        try:
            si = int(start) if start is not None else None
            ei = int(end)   if end   is not None else si
        except (TypeError, ValueError):
            continue
        if si is None:
            continue

        if ftype == "active site" or "catalytic" in text:
            key = ("Catalytic site", si, ei)
            if key in seen_s:
                continue
            seen_s.add(key)
            sites.append({"label": "Catalytic site", "start": si, "end": ei,
                          "description": desc or "Catalytic residue",
                          "color": _FUNC_COLORS["Catalytic site"]})
            continue

        if ftype not in {"domain", "region"}:
            continue

        if "exonuclease" in text:
            lbl = "Exonuclease"
        elif "polymerase" in text:
            lbl = "Polymerase core"
        else:
            continue

        key = (lbl, si, ei)
        if key in seen_r:
            continue
        seen_r.add(key)
        regions.append({"label": desc or lbl, "start": si, "end": ei,
                        "description": desc or lbl,
                        "color": _FUNC_COLORS[lbl]})

    return regions, sites


def _residues_in_range(
    pdb_coords: dict[int, dict], start: int, end: int
) -> list[tuple[int, dict]]:
    lo, hi = sorted((start, end))
    return [(r, pdb_coords[r]) for r in sorted(pdb_coords) if lo <= r <= hi]


def _centroid(items: list[tuple[int, dict]]) -> tuple[float, float, float]:
    if not items:
        return 0.0, 0.0, 0.0
    n = float(len(items))
    return (
        sum(c["x"] for _, c in items) / n,
        sum(c["y"] for _, c in items) / n,
        sum(c["z"] for _, c in items) / n,
    )


def _add_functional_overlays(
    fig: go.Figure,
    pdb_coords: dict[int, dict],
    annotations: list[dict] | None,
) -> None:
    regions, sites = _extract_functional_annotations(annotations)

    for region in regions:
        pts = _residues_in_range(pdb_coords, region["start"], region["end"])
        if len(pts) < 2:
            continue
        fig.add_trace(go.Scatter3d(
            x=[c["x"] for _, c in pts],
            y=[c["y"] for _, c in pts],
            z=[c["z"] for _, c in pts],
            mode="lines",
            line={"color": region["color"], "width": 10},
            name=region["label"], showlegend=False,
            hovertemplate=(
                f"<b>{region['label']}</b><br>"
                f"Residues: {region['start']}-{region['end']}<br>"
                f"{region['description']}<extra></extra>"
            ),
        ))
        cx, cy, cz = _centroid(pts)
        fig.add_trace(go.Scatter3d(
            x=[cx], y=[cy], z=[cz],
            mode="text", text=[region["label"]],
            textposition="top center",
            textfont={"color": region["color"], "size": 12},
            showlegend=False, hoverinfo="skip",
        ))

    if not sites:
        return

    site_pts = [
        (s, _centroid(_residues_in_range(pdb_coords, s["start"], s["end"])))
        for s in sites
        if _residues_in_range(pdb_coords, s["start"], s["end"])
    ]
    if not site_pts:
        return

    fig.add_trace(go.Scatter3d(
        x=[p[1][0] for p in site_pts],
        y=[p[1][1] for p in site_pts],
        z=[p[1][2] for p in site_pts],
        mode="markers",
        marker={"size": 9, "symbol": "cross",
                "color": _FUNC_COLORS["Catalytic site"],
                "line": {"color": "#ffffff", "width": 1.5}},
        text=[
            f"{s['description']} ({s['start']}"
            + ("" if s["start"] == s["end"] else f"-{s['end']}") + ")"
            for s, _ in site_pts
        ],
        name="Catalytic site", showlegend=False,
        hovertemplate="<b>%{text}</b><extra></extra>",
    ))


# ---------------------------------------------------------------------------
# Public: build_3d_fingerprint
# ---------------------------------------------------------------------------

_GEN_COLORS: list[str] = [
    "#a78bfa", "#fb923c", "#22d3ee", "#4ade80",
    "#f472b6", "#facc15", "#34d399", "#60a5fa",
    "#fb7185", "#e879f9",
]


def build_3d_fingerprint(
    lineage_mutations: list[dict[str, Any]],
    backbone: list[dict[str, Any]],
    pdb_coords: dict[int, dict[str, Any]],
    variant_id: int | float,
    uniprot_id: str | None = None,
    structure_source: str = "Structure",
    feature_annotations: list[dict[str, Any]] | None = None,
) -> go.Figure:
    """
    Build a structure-first 3D mutation map.

    Backbone is pLDDT-coloured (red=low confidence → blue=high).
    Mutations are plotted as generation-coloured spheres
    (diamond = non-synonymous, circle = synonymous).
    Functional region overlays are added when UniProt annotations are supplied.

    Parameters
    ----------
    lineage_mutations   list of mutation dicts — each must have:
                        position, wt_aa, mut_aa, mutation_type, generation,
                        aa_change
    backbone            list[{r, x, y, z, plddt}]
    pdb_coords          dict[int, {x, y, z, plddt}]
    variant_id          numeric variant index (title)
    uniprot_id          UniProt accession (subtitle)
    structure_source    source label (subtitle)
    feature_annotations UniProt feature list for overlay (optional)
    """
    fig = go.Figure()

    # Backbone ribbon
    if backbone:
        fig.add_trace(go.Scatter3d(
            x=[r["x"] for r in backbone],
            y=[r["y"] for r in backbone],
            z=[r["z"] for r in backbone],
            mode="lines",
            line={
                "color": [r.get("plddt", 70.0) for r in backbone],
                "colorscale": [
                    [0.00, "#f87171"],
                    [0.35, "#facc15"],
                    [0.55, "#38bdf8"],
                    [1.00, "#2563eb"],
                ],
                "cmin": 0, "cmax": 100, "width": 5,
            },
            hoverinfo="skip", showlegend=False, name="C-alpha backbone",
        ))
        _add_functional_overlays(fig, pdb_coords, feature_annotations)

    generations = sorted({m["generation"] for m in lineage_mutations})

    _STYLES = (
        ("non-synonymous", "Non-synonymous", 9,  "diamond"),
        ("synonymous",     "Synonymous",     6,  "circle"),
    )

    for gen in generations:
        gen_color = _GEN_COLORS[(gen - 1) % len(_GEN_COLORS)]
        for mut_key, mut_label, m_size, m_symbol in _STYLES:
            xs: list[float] = []
            ys: list[float] = []
            zs: list[float] = []
            custom: list = []

            for m in lineage_mutations:
                if m["generation"] != gen:
                    continue
                t = ("synonymous"
                     if m["mutation_type"] == "synonymous"
                     else "non-synonymous")
                if t != mut_key:
                    continue
                coord = pdb_coords.get(int(m["position"]))
                if coord is None:
                    continue
                xs.append(coord["x"])
                ys.append(coord["y"])
                zs.append(coord["z"])
                custom.append([
                    m.get("aa_change",
                          f"{m['wt_aa']}{m['position']}{m['mut_aa']}"),
                    m["position"],
                    m["wt_aa"],
                    m["mut_aa"],
                    m["mutation_type"],
                    gen,
                    round(coord.get("plddt", 70.0), 1),
                ])

            if not xs:
                continue

            fig.add_trace(go.Scatter3d(
                x=xs, y=ys, z=zs,
                mode="markers",
                name=f"Gen {gen} — {mut_label}",
                showlegend=True,
                meta={"generation": gen, "mutation_type": mut_key},
                marker={
                    "symbol": m_symbol, "size": m_size,
                    "color": gen_color, "opacity": 0.90,
                    "line": {"color": "rgba(0,0,0,0.30)", "width": 0.5},
                },
                customdata=custom,
                hovertemplate=(
                    "<b>%{customdata[0]}</b><br>"
                    "Position: %{customdata[1]}<br>"
                    "Change: %{customdata[2]} → %{customdata[3]}<br>"
                    "Type: %{customdata[4]}<br>"
                    "Generation: %{customdata[5]}<br>"
                    "pLDDT: %{customdata[6]}<extra></extra>"
                ),
            ))

    subtitle = (
        f"UniProt {uniprot_id} · {structure_source}"
        if uniprot_id else structure_source
    )

    fig.update_layout(
        title={
            "text": (
                f"Variant {variant_id} — 3D Mutation Map<br>"
                f"<sup>{subtitle}</sup>"
            ),
            "font": {"size": 15},
        },
        scene={
            "xaxis": {"visible": False, "showgrid": False, "zeroline": False},
            "yaxis": {"visible": False, "showgrid": False, "zeroline": False},
            "zaxis": {"visible": False, "showgrid": False, "zeroline": False},
            "bgcolor": "rgba(6,9,15,1)",
            "camera": {"eye": {"x": 1.5, "y": 1.5, "z": 0.8}},
            "aspectmode": "data",
        },
        paper_bgcolor="rgba(6,9,15,1)",
        plot_bgcolor="rgba(6,9,15,1)",
        font={"color": "#d1dce8"},
        legend={
            "title": {"text": "Generation / Type", "font": {"size": 11}},
            "bgcolor": "rgba(17,24,39,0.9)",
            "bordercolor": "#1e3050", "borderwidth": 1,
            "itemclick": "toggle",
            "itemdoubleclick": "toggleothers",
        },
        margin={"l": 0, "r": 0, "t": 80, "b": 0},
        height=760,
        hoverlabel={
            "bgcolor": "#162032", "bordercolor": "#1e3050",
            "font": {"family": "monospace", "size": 12, "color": "#d1dce8"},
        },
    )

    return fig
