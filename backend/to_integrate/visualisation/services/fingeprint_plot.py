"""
fingerprint_plot.py
-------------------
Part 4a — Visualisation Layer
================================
Builds the circular mutation fingerprint Plotly figure and, optionally,
maps mutated residues onto a 3D protein structure (AlphaFill / PDB).

This module has ONE job: take an *analysis* dict produced by
mutation_analysis.analyze_target_variant() and return a Plotly Figure.
It never touches sequence data directly.

Public API
----------
    resolve_structure(uniprot_id, pdb_file, no_3d) -> (coords, info)
    build_plotly_figure(analysis, structure_coords, structure_info) -> go.Figure

CLI usage
---------
    python fingerprint_plot.py --variant-id 295 --no-3d
    python fingerprint_plot.py --variant-id 295 --uniprot-id O34996
    python fingerprint_plot.py --variant-id 295 \\
        --input DE_BSU_Pol_Batch_1.json \\
        --mutation-csv variant_295_mutations.csv

Dependencies
------------
    plotly        (pip install plotly)
    requests      (pip install requests)   [only for 3D structure download]
    mutation_analysis.py  (same directory)
"""

from __future__ import annotations

import argparse
import logging
import math
import tempfile
from pathlib import Path
from typing import Any

import requests
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Local analysis module — same package
from .mutation_analysis import (
    read_fasta,
    load_variants_table,
    analyze_target_variant,
    get_variant_lineage,
    analyze_lineage_mutations,
    write_mutation_csv,
    infer_uniprot_id,
)

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# 3-D structure helpers
# ---------------------------------------------------------------------------

def _fetch_url_text(url: str, timeout: float = 20.0) -> str:
    resp = requests.get(url, timeout=timeout)
    resp.raise_for_status()
    return resp.text


def _fetch_alphafill_model(uniprot_id: str, cache_dir: Path) -> tuple[Path, str]:
    """Try AlphaFill versions v4 → v1 and return the first model text that works."""
    cache_dir.mkdir(parents=True, exist_ok=True)

    for version in ("v4", "v3", "v2", "v1"):
        model_id = f"AF-{uniprot_id}-F1-model_{version}"
        filename = f"{model_id}.cif"
        cached = cache_dir / filename
        if cached.exists():
            return cached, f"AlphaFill ({model_id})"

        url = f"https://alphafill.eu/v1/aff/{model_id}"
        try:
            text = _fetch_url_text(url)
        except Exception:
            continue

        if "_atom_site." not in text:
            continue

        cached.write_text(text, encoding="utf-8")
        return cached, f"AlphaFill ({model_id})"

    raise FileNotFoundError(f"No AlphaFill model found for UniProt ID '{uniprot_id}'")


def _fetch_uniprot_primary_pdb(uniprot_id: str, cache_dir: Path) -> tuple[Path, str]:
    """Fetch the first available PDB structure cross-referenced in UniProt."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    data = requests.get(
        f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json", timeout=20.0
    ).json()

    pdb_ids = [
        str(ref["id"])
        for ref in data.get("uniProtKBCrossReferences", [])
        if ref.get("database") == "PDB" and ref.get("id")
    ]

    if not pdb_ids:
        raise FileNotFoundError(f"No PDB cross-references in UniProt entry '{uniprot_id}'")

    for pdb_id in pdb_ids:
        cache_path = cache_dir / f"{pdb_id}.pdb"
        if cache_path.exists():
            return cache_path, f"PDB ({pdb_id})"
        try:
            text = _fetch_url_text(f"https://files.rcsb.org/download/{pdb_id}.pdb")
        except Exception:
            continue
        if "ATOM" not in text:
            continue
        cache_path.write_text(text, encoding="utf-8")
        return cache_path, f"PDB ({pdb_id})"

    raise FileNotFoundError(f"Could not download any PDB structure for UniProt ID '{uniprot_id}'")


def _parse_pdb_ca_coordinates(
    pdb_path: Path,
) -> tuple[dict[int, tuple[float, float, float]], str]:
    """
    Parse Cα atom coordinates from a PDB file.

    Returns:
        coords  {residue_number -> (x, y, z)}
        chain   chain ID of the largest chain
    """
    chains: dict[str, dict[int, tuple[float, float, float]]] = {}

    for line in pdb_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith("ATOM"):
            continue
        if line[12:16].strip() != "CA":
            continue
        chain_id = line[21].strip() or "_"
        try:
            pos = int(line[22:26].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
        except Exception:
            continue
        chains.setdefault(chain_id, {})[pos] = (x, y, z)

    if not chains:
        raise ValueError(f"No Cα atoms parsed from {pdb_path}")

    best_chain, best_coords = max(chains.items(), key=lambda kv: len(kv[1]))
    return best_coords, best_chain


def _parse_mmcif_ca_coordinates(
    structure_path: Path,
) -> tuple[dict[int, tuple[float, float, float]], str]:
    """
    Parse C-alpha coordinates from an mmCIF/ModelCIF atom_site loop.

    AlphaFill serves structures in ModelCIF text, not PDB, so the parser
    extracts the relevant atom_site columns by name rather than relying on
    fixed-width PDB columns.
    """
    lines = structure_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    atom_columns: list[str] = []
    data_start: int | None = None

    for i, line in enumerate(lines):
        if line.strip() != "loop_":
            continue
        columns: list[str] = []
        j = i + 1
        while j < len(lines) and lines[j].lstrip().startswith("_"):
            columns.append(lines[j].strip())
            j += 1
        if columns and columns[0].startswith("_atom_site."):
            atom_columns = columns
            data_start = j
            break

    if not atom_columns or data_start is None:
        raise ValueError(f"No atom_site loop parsed from {structure_path}")

    col_index = {name: idx for idx, name in enumerate(atom_columns)}

    def _pick_index(*names: str) -> int:
        for name in names:
            if name in col_index:
                return col_index[name]
        raise KeyError(names[0])

    atom_idx = _pick_index("_atom_site.auth_atom_id", "_atom_site.label_atom_id")
    chain_idx = _pick_index("_atom_site.auth_asym_id", "_atom_site.label_asym_id")
    seq_idx = _pick_index("_atom_site.auth_seq_id", "_atom_site.label_seq_id")
    x_idx = _pick_index("_atom_site.Cartn_x")
    y_idx = _pick_index("_atom_site.Cartn_y")
    z_idx = _pick_index("_atom_site.Cartn_z")

    chains: dict[str, dict[int, tuple[float, float, float]]] = {}

    for line in lines[data_start:]:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped == "#":
            if chains:
                break
            continue
        if stripped == "loop_" or stripped.startswith("_"):
            if chains:
                break
            continue

        parts = stripped.split()
        if len(parts) < len(atom_columns):
            continue

        try:
            atom_name = parts[atom_idx]
            if atom_name != "CA":
                continue
            chain_id = parts[chain_idx]
            if chain_id in {".", "?"}:
                chain_id = "_"
            seq_text = parts[seq_idx]
            if seq_text in {".", "?"}:
                continue
            pos = int(seq_text)
            x = float(parts[x_idx])
            y = float(parts[y_idx])
            z = float(parts[z_idx])
        except Exception:
            continue

        chains.setdefault(chain_id, {})[pos] = (x, y, z)

    if not chains:
        raise ValueError(f"No C-alpha atoms parsed from {structure_path}")

    best_chain, best_coords = max(chains.items(), key=lambda kv: len(kv[1]))
    return best_coords, best_chain


def _parse_structure_ca_coordinates(
    structure_path: Path,
) -> tuple[dict[int, tuple[float, float, float]], str]:
    text = structure_path.read_text(encoding="utf-8", errors="ignore")
    if "_atom_site." in text and not text.lstrip().startswith("ATOM"):
        return _parse_mmcif_ca_coordinates(structure_path)
    return _parse_pdb_ca_coordinates(structure_path)


def resolve_structure(
    uniprot_id: str | None,
    pdb_file: Path | None,
    no_3d: bool,
) -> tuple[dict[int, tuple[float, float, float]] | None, dict[str, str]]:
    """
    Resolve Cα coordinates for 3D mapping.

    Priority order:
      1. --no-3d flag → skip entirely
      2. --pdb-file path → load local PDB
      3. AlphaFill download (requires uniprot_id)
      4. UniProt PDB cross-reference download
      5. Graceful fallback with status message

    Returns:
        (coords_dict, info_dict)
        coords_dict is None when 3D mapping is unavailable or disabled.
    """
    info: dict[str, str] = {}

    if no_3d:
        info["status"] = "3D mapping disabled via --no-3d"
        return None, info

    if pdb_file is not None:
        if not pdb_file.exists():
            raise FileNotFoundError(f"Provided --pdb-file not found: {pdb_file}")
        coords, chain = _parse_structure_ca_coordinates(pdb_file)
        info.update(status="3D mapping from local PDB", source=str(pdb_file), chain=chain)
        return coords, info

    if not uniprot_id:
        info["status"] = "No UniProt ID available — skipping 3D mapping"
        return None, info

    cache_dir = Path(tempfile.gettempdir()) / "fingerprint_plot_structures"

    # Try AlphaFill first
    try:
        pdb_path, label = _fetch_alphafill_model(uniprot_id, cache_dir)
        coords, chain = _parse_structure_ca_coordinates(pdb_path)
        info.update(status="3D mapping enabled", source=label, chain=chain)
        return coords, info
    except Exception as exc:
        logger.warning("AlphaFill retrieval failed for %s: %s", uniprot_id, exc)

    # Fallback to UniProt PDB cross-reference
    try:
        pdb_path, label = _fetch_uniprot_primary_pdb(uniprot_id, cache_dir)
        coords, chain = _parse_structure_ca_coordinates(pdb_path)
        info.update(status="3D mapping enabled", source=label, chain=chain)
        return coords, info
    except Exception as exc:
        logger.warning("UniProt/PDB retrieval failed for %s: %s", uniprot_id, exc)

    info["status"] = "3D mapping unavailable (download or parsing failed)"
    return None, info


# ---------------------------------------------------------------------------
# Plotly figure builder
# ---------------------------------------------------------------------------

def build_plotly_figure(
    analysis: dict[str, Any],
    structure_coords: dict[int, tuple[float, float, float]] | None,
    structure_info: dict[str, str],
) -> go.Figure:
    """
    Build the circular mutation fingerprint Plotly figure.

    Layout
    ------
    * Left panel (always present):
        - Outer grey circle  = WT protein backbone
        - Inner blue circle  = Variant protein backbone
        - Radial connectors  = Mutation positions
        - Red  triangles     = Non-synonymous mutations
        - Blue circles       = Synonymous mutations
        - Residue tick labels every 100 aa (or 250 aa for long proteins)
    * Right panel (when 3D structure is available):
        - Cα backbone trace of the resolved structure
        - Colour-coded mutation dots mapped onto structure residues

    Parameters
    ----------
    analysis        dict returned by mutation_analysis.analyze_target_variant()
    structure_coords  {residue_number -> (x, y, z)} or None
    structure_info    metadata dict from resolve_structure()

    Returns
    -------
    plotly.graph_objects.Figure
    """
    mutations = analysis["mutations"]
    alignment_length = max(int(analysis["alignment_length"]), 1)
    wt_len = max(len(analysis["wt_protein_sequence"]), 1)

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def angle_from_aligned_pos(aligned_pos: int) -> float:
        """Convert an alignment column to radians (12-o'clock start, clockwise)."""
        return (2.0 * math.pi * (int(aligned_pos) - 1) / alignment_length) - (math.pi / 2.0)

    def circle_trace(radius: float, color: str, width: int, name: str) -> go.Scatter:
        xs = [radius * math.cos(math.radians(d)) for d in range(361)]
        ys = [radius * math.sin(math.radians(d)) for d in range(361)]
        return go.Scatter(
            x=xs, y=ys,
            mode="lines",
            line={"color": color, "width": width},
            hoverinfo="skip",
            name=name,
        )

    # ------------------------------------------------------------------
    # Figure skeleton
    # ------------------------------------------------------------------

    has_3d = structure_coords is not None
    if has_3d:
        fig = make_subplots(
            rows=1, cols=2,
            specs=[[{"type": "xy"}, {"type": "scene"}]],
            column_widths=[0.56, 0.44],
            subplot_titles=(
                "Circular WT vs Variant Mutation Fingerprint",
                f"3D Structure Mutation Map ({structure_info.get('source', 'structure')})",
            ),
            horizontal_spacing=0.05,
        )
    else:
        fig = go.Figure()

    # ------------------------------------------------------------------
    # 2D circular traces
    # ------------------------------------------------------------------

    WT_RADIUS = 1.0
    VAR_RADIUS = 0.78

    traces_2d: list[go.BaseTraceType] = [
        circle_trace(WT_RADIUS, "#94a3b8", 20, "WT circular protein"),
        circle_trace(VAR_RADIUS, "#3b82f6", 16, "Variant mapped protein"),
    ]

    # Connector lines from WT ring to variant ring at each mutation
    conn_x: list[float] = []
    conn_y: list[float] = []

    # Separate lists for synonymous and non-synonymous markers
    syn_x: list[float] = []
    syn_y: list[float] = []
    syn_custom: list[list[Any]] = []

    nonsyn_x: list[float] = []
    nonsyn_y: list[float] = []
    nonsyn_custom: list[list[Any]] = []

    for m in mutations:
        aln_pos = int(m.get("aligned_position", m["position"]))
        theta = angle_from_aligned_pos(aln_pos)

        wt_x = WT_RADIUS * math.cos(theta)
        wt_y = WT_RADIUS * math.sin(theta)
        var_x = VAR_RADIUS * math.cos(theta)
        var_y = VAR_RADIUS * math.sin(theta)

        conn_x.extend([wt_x, var_x, None])
        conn_y.extend([wt_y, var_y, None])

        row_custom = [
            m["position"],
            m.get("aligned_position", m["position"]),
            m["wt_aa"],
            m["mut_aa"],
            m["mutation_type"],
            m["wt_codon"],
            m["mut_codon"],
            m["aa_change"],
        ]

        if m["mutation_type"] == "synonymous":
            syn_x.append(var_x)
            syn_y.append(var_y)
            syn_custom.append(row_custom)
        else:
            nonsyn_x.append(var_x)
            nonsyn_y.append(var_y)
            nonsyn_custom.append(row_custom)

    hover_template = (
        "WT position %{customdata[0]} (aligned %{customdata[1]})<br>"
        "AA change: %{customdata[2]}→%{customdata[3]} (%{customdata[4]})<br>"
        "Codon: %{customdata[5]}→%{customdata[6]}<br>"
        "Label: %{customdata[7]}<extra></extra>"
    )

    traces_2d.append(go.Scatter(
        x=conn_x, y=conn_y,
        mode="lines",
        line={"color": "rgba(71,85,105,0.30)", "width": 1},
        hoverinfo="skip",
        showlegend=False,
        name="WT-to-variant mapping",
    ))

    traces_2d.append(go.Scatter(
        x=nonsyn_x, y=nonsyn_y,
        mode="markers",
        name="Non-synonymous",
        marker={
            "size": 12,
            "color": "#dc2626",
            "symbol": "triangle-down",
            "line": {"color": "#1f2937", "width": 1},
        },
        customdata=nonsyn_custom,
        hovertemplate=hover_template,
    ))

    traces_2d.append(go.Scatter(
        x=syn_x, y=syn_y,
        mode="markers",
        name="Synonymous",
        marker={
            "size": 10,
            "color": "#0ea5e9",
            "symbol": "circle",
            "line": {"color": "#1f2937", "width": 1},
        },
        customdata=syn_custom,
        hovertemplate=hover_template,
    ))

    # Residue-position tick labels on the outer ring
    tick_step = 250 if wt_len >= 1000 else 100
    tick_x: list[float] = []
    tick_y: list[float] = []
    tick_text: list[str] = []
    for pos in range(1, wt_len + 1, tick_step):
        theta = (2.0 * math.pi * (pos - 1) / wt_len) - (math.pi / 2.0)
        tick_x.append(1.16 * math.cos(theta))
        tick_y.append(1.16 * math.sin(theta))
        tick_text.append(str(pos))

    traces_2d.append(go.Scatter(
        x=tick_x, y=tick_y,
        mode="text",
        text=tick_text,
        textfont={"size": 10, "color": "#334155"},
        hoverinfo="skip",
        showlegend=False,
        name="Position labels",
    ))

    # Add 2D traces to the figure
    if has_3d:
        for tr in traces_2d:
            fig.add_trace(tr, row=1, col=1)
    else:
        for tr in traces_2d:
            fig.add_trace(tr)

    # ------------------------------------------------------------------
    # 3D structure traces (optional)
    # ------------------------------------------------------------------

    mapped_count = 0
    missing_positions: list[int] = []

    if has_3d:
        assert structure_coords is not None

        residue_positions = sorted(structure_coords.keys())
        backbone_x = [structure_coords[p][0] for p in residue_positions]
        backbone_y = [structure_coords[p][1] for p in residue_positions]
        backbone_z = [structure_coords[p][2] for p in residue_positions]

        fig.add_trace(
            go.Scatter3d(
                x=backbone_x, y=backbone_y, z=backbone_z,
                mode="lines",
                line={"color": "#94a3b8", "width": 6},
                name="Structure backbone (Cα)",
                hoverinfo="skip",
                showlegend=False,
            ),
            row=1, col=2,
        )

        mut3d_x: list[float] = []
        mut3d_y: list[float] = []
        mut3d_z: list[float] = []
        mut3d_custom: list[list[Any]] = []
        mut3d_colors: list[str] = []

        for m in mutations:
            pos = int(m["position"])
            if pos not in structure_coords:
                missing_positions.append(pos)
                continue
            x, y, z = structure_coords[pos]
            mapped_count += 1
            mut3d_x.append(x)
            mut3d_y.append(y)
            mut3d_z.append(z)
            mut3d_custom.append([
                m["position"], m["wt_aa"], m["mut_aa"],
                m["mutation_type"], m["aa_change"],
            ])
            mut3d_colors.append("#dc2626" if m["mutation_type"] != "synonymous" else "#0ea5e9")

        if mut3d_x:
            fig.add_trace(
                go.Scatter3d(
                    x=mut3d_x, y=mut3d_y, z=mut3d_z,
                    mode="markers",
                    name="Mapped variant mutations",
                    marker={
                        "size": 6,
                        "color": mut3d_colors,
                        "line": {"color": "#111827", "width": 0.5},
                    },
                    customdata=mut3d_custom,
                    hovertemplate=(
                        "Position %{customdata[0]}<br>"
                        "AA %{customdata[1]}→%{customdata[2]} (%{customdata[3]})<br>"
                        "Change %{customdata[4]}<extra></extra>"
                    ),
                ),
                row=1, col=2,
            )

    # ------------------------------------------------------------------
    # Layout
    # ------------------------------------------------------------------

    title = (
        f"Variant {analysis['variant_id']} Mutation Fingerprint "
        f"(n={analysis['num_mutations']} mutations vs WT)"
    )

    fig.update_layout(
        title=title,
        template="plotly_white",
        height=860,
        margin={"l": 30, "r": 30, "t": 88, "b": 24},
        legend={"orientation": "h", "x": 0.02, "y": 1.02},
    )

    if has_3d:
        fig.update_xaxes(visible=False, range=[-1.30, 1.30], row=1, col=1)
        fig.update_yaxes(
            visible=False, range=[-1.30, 1.30],
            scaleanchor="x", scaleratio=1, row=1, col=1,
        )
        fig.update_scenes(
            xaxis_visible=False, yaxis_visible=False, zaxis_visible=False,
            bgcolor="rgba(0,0,0,0)", aspectmode="data",
            row=1, col=2,
        )
    else:
        fig.update_xaxes(visible=False, range=[-1.30, 1.30])
        fig.update_yaxes(visible=False, range=[-1.30, 1.30], scaleanchor="x", scaleratio=1)

    # Summary annotation (bottom-left)
    syn = analysis["num_synonymous"]
    nonsyn = analysis["num_nonsynonymous"]
    fig.add_annotation(
        xref="paper", yref="paper",
        x=0.01, y=0.01,
        showarrow=False, align="left",
        bordercolor="#cbd5e1", borderwidth=1, bgcolor="#f8fafc",
        font={"size": 12, "color": "#1f2937"},
        text=(
            f"Variant {analysis['variant_id']} | "
            f"Generation {analysis['generation']} | "
            f"Parent {analysis['parent_variant_id']}<br>"
            f"WT length: {len(analysis['wt_protein_sequence'])} aa | "
            f"Variant length: {len(analysis['variant_protein_sequence'])} aa<br>"
            f"Mutations: {analysis['num_mutations']} "
            f"(non-synonymous: {nonsyn}, synonymous: {syn})<br>"
            f"Alignment — matches: {analysis['alignment_stats']['matches']}, "
            f"mismatches: {analysis['alignment_stats']['mismatches']}, "
            f"gaps: {analysis['alignment_stats']['gaps']}"
        ),
    )

    # 3D annotation (bottom-right)
    if has_3d:
        preview = ", ".join(str(p) for p in missing_positions[:12])
        more = " ..." if len(missing_positions) > 12 else ""
        fig.add_annotation(
            xref="paper", yref="paper",
            x=0.99, y=0.01,
            showarrow=False, align="right",
            bordercolor="#cbd5e1", borderwidth=1, bgcolor="#f8fafc",
            font={"size": 12, "color": "#1f2937"},
            text=(
                f"3D source: {structure_info.get('source', 'n/a')} | "
                f"chain {structure_info.get('chain', '?')}<br>"
                f"Mapped mutations in 3D: {mapped_count}/{analysis['num_mutations']}<br>"
                f"Missing residues in structure: {len(missing_positions)}"
                + (f" ({preview}{more})" if missing_positions else "")
            ),
        )
    else:
        fig.add_annotation(
            xref="paper", yref="paper",
            x=0.99, y=0.01,
            showarrow=False, align="right",
            bordercolor="#cbd5e1", borderwidth=1, bgcolor="#f8fafc",
            font={"size": 12, "color": "#1f2937"},
            text=structure_info.get("status", "3D mapping unavailable"),
        )

    return fig


# ---------------------------------------------------------------------------
# Linear (unfolded) mutation fingerprint  — matches Figure 3 in the brief
# ---------------------------------------------------------------------------

# One colour per generation (up to 10).  Matches the Figure 3 legend order.
_GENERATION_COLORS: list[str] = [
    "#f97316",  # Gen 1  — orange
    "#22c55e",  # Gen 2  — green
    "#ef4444",  # Gen 3  — red
    "#a855f7",  # Gen 4  — purple
    "#8b5cf6",  # Gen 5  — violet
    "#64748b",  # Gen 6  — slate
    "#84cc16",  # Gen 7  — lime
    "#d4d4aa",  # Gen 8  — tan
    "#a3e635",  # Gen 9  — yellow-green
    "#3b82f6",  # Gen 10 — blue
]


def _assign_label_rows(positions: list[int], min_gap: int = 30) -> list[int]:
    """
    Assign a vertical row index (0 = lowest) to each mutation so that labels
    at nearby positions do not collide.

    Uses a greedy left-to-right sweep: each mutation is placed in the lowest
    row where the previous occupant is at least *min_gap* residues to the
    left.
    """
    row_last: list[int] = []  # rightmost position currently in each row
    rows: list[int] = []
    for pos in positions:
        placed = False
        for idx, last in enumerate(row_last):
            if pos - last >= min_gap:
                row_last[idx] = pos
                rows.append(idx)
                placed = True
                break
        if not placed:
            rows.append(len(row_last))
            row_last.append(pos)
    return rows


def build_linear_fingerprint(
    lineage_mutations: list[dict[str, Any]],
    wt_protein_len: int,
    variant_id: int,
) -> go.Figure:
    """
    Build a linear (unfolded) mutation fingerprint Plotly figure.

    Layout (mirrors Figure 3 in the brief)
    ----------------------------------------
    * Horizontal grey bar  = WT protein backbone
    * Downward triangles   = mutation positions, one per AA change
    * Colour               = generation in which the mutation was introduced
    * Text labels above    = 'WT>Mut (m)' or '(s)' for non/synonymous
    * Labels are stacked vertically at nearby positions to avoid overlap
    * X-axis               = Amino Acid Position (1-based)

    Parameters
    ----------
    lineage_mutations   list returned by analyze_lineage_mutations()
    wt_protein_len      length of the WT protein in residues
    variant_id          used only for the figure title
    """
    # Show inline text labels only when mutation count is manageable.
    # Above the threshold the labels collapse into an unreadable wall and
    # push the backbone off-screen; hover text carries the same information.
    _LABEL_THRESHOLD = 80
    show_labels = len(lineage_mutations) <= _LABEL_THRESHOLD

    # Vertical layout — mirrors Figure 3:
    #
    #   label text  (stacked rows above the triangles)
    #       |        (thin connector line, drawn behind the triangle)
    #       ▼        (triangle marker, pointing DOWN toward backbone)
    #   ━━━━━━━━━   (grey protein backbone bar)
    #   (empty)     (breathing room below)
    #
    ROW_HEIGHT   = 0.55   # gap between label rows
    MARKER_Y     = 0.35   # triangles sit just above the backbone
    LABEL_Y_BASE = 0.75   # y of the lowest (row-0) label centre
    BACKBONE_TOP = 0.12   # top edge of the backbone rect
    Y_BOTTOM     = -0.80  # empty space below backbone

    fig = go.Figure()

    # --- Protein backbone bar ---
    fig.add_shape(
        type="rect",
        x0=0, x1=wt_protein_len,
        y0=-0.12, y1=BACKBONE_TOP,
        fillcolor="#94a3b8",
        line={"width": 0},
        layer="below",
    )

    # --- One scatter trace per generation (for clean legend) ---
    generations = sorted({m["generation"] for m in lineage_mutations})

    # Pre-assign label rows only when we will actually show them.
    if show_labels:
        all_sorted = sorted(lineage_mutations, key=lambda m: m["position"])
        row_indices = _assign_label_rows([m["position"] for m in all_sorted])
        row_lookup: dict[tuple[int, int], int] = {
            (m["generation"], m["position"]): r
            for m, r in zip(all_sorted, row_indices)
        }
        max_row = max(row_indices) if row_indices else 0
    else:
        row_indices = []
        row_lookup = {}
        max_row = 0

    for gen in generations:
        gen_muts = [m for m in lineage_mutations if m["generation"] == gen]
        color = _GENERATION_COLORS[(gen - 1) % len(_GENERATION_COLORS)]

        xs: list[float] = []
        ys_text: list[float] = []
        labels: list[str] = []
        hover_texts: list[str] = []

        for m in gen_muts:
            pos = int(m["position"])
            mut_type = m["mutation_type"]
            if mut_type == "insertion":
                suffix = "(i)"
            elif mut_type == "deletion":
                suffix = "(d)"
            elif mut_type == "synonymous":
                suffix = "(s)"
            else:
                suffix = "(m)"
            label = f"{m['wt_aa']}>{m['mut_aa']}<br>{suffix}"
            row = row_lookup.get((gen, pos), 0)

            xs.append(pos)
            ys_text.append(LABEL_Y_BASE + row * ROW_HEIGHT)
            labels.append(label)
            codon_line = (
                f"Codon: {m['wt_codon']}→{m['mut_codon']}<br>"
                if m.get("wt_codon", "n/a") != "n/a" else ""
            )
            hover_texts.append(
                f"<b>{m['aa_change']}</b><br>"
                f"Position: {pos}<br>"
                f"Type: {m['mutation_type']}<br>"
                + codon_line
                + f"Generation: {gen}"
            )

        # Triangles — always shown, sitting just above the backbone
        fig.add_trace(go.Scatter(
            x=xs, y=[MARKER_Y] * len(xs),
            mode="markers",
            name=f"Generation {gen}",
            marker={
                "symbol": "triangle-down",
                "size": 10 if show_labels else 7,
                "color": color,
                "line": {"color": "#1f2937", "width": 1},
            },
            hovertext=hover_texts,
            hovertemplate="%{hovertext}<extra></extra>",
            showlegend=True,
        ))

        if show_labels:
            # Text labels above the triangles
            fig.add_trace(go.Scatter(
                x=xs, y=ys_text,
                mode="text",
                text=labels,
                textfont={"size": 9, "color": color},
                hoverinfo="skip",
                showlegend=False,
                name=f"Generation {gen} labels",
            ))
            # Connector line: runs from backbone top up through triangle to label
            for x, yt in zip(xs, ys_text):
                fig.add_shape(
                    type="line",
                    x0=x, x1=x,
                    y0=BACKBONE_TOP,
                    y1=yt - 0.12,       # stop just below the label text
                    line={"color": "#94a3b8", "width": 1},
                    layer="below",      # drawn behind the triangle marker
                )

    # --- Axis / layout ---
    if show_labels:
        y_top = LABEL_Y_BASE + (max_row + 1) * ROW_HEIGHT + 0.2
        fig_h = max(480, 300 + max_row * 60)
        note  = ""
    else:
        y_top = MARKER_Y + 1.2   # compact: a little above the triangles
        fig_h = 420
        note  = f"  —  {len(lineage_mutations)} mutations (hover a triangle for details)"

    fig.update_layout(
        title=f"Variant {variant_id} — Linear Mutation Fingerprint by Generation{note}",
        xaxis={
            "title": "Amino Acid Position",
            "range": [-15, wt_protein_len + 15],
            "showgrid": True,
            "gridcolor": "#e2e8f0",
            "tickmode": "linear",
            "dtick": 100,
            "zeroline": False,
        },
        yaxis={
            "visible": False,
            "range": [Y_BOTTOM, y_top],
        },
        template="plotly_white",
        height=fig_h,
        legend={
            "title": {"text": "Generation"},
            "orientation": "v",
            "x": 1.01, "y": 1,
            "bgcolor": "rgba(255,255,255,0.9)",
            "bordercolor": "#cbd5e1", "borderwidth": 1,
        },
        margin={"l": 50, "r": 160, "t": 60, "b": 60},
        hovermode="closest",
    )

    return fig


# ---------------------------------------------------------------------------
# 3-D mutation fingerprint  — mutation-type coloured spheres
# ---------------------------------------------------------------------------

_GEN_COLORS_3D: list[str] = [
    "#a78bfa",
    "#fb923c",
    "#22d3ee",
    "#4ade80",
    "#f472b6",
    "#facc15",
    "#34d399",
    "#60a5fa",
    "#fb7185",
    "#e879f9",
]


_FUNCTIONAL_REGION_COLORS: dict[str, str] = {
    "Exonuclease": "#14b8a6",
    "Polymerase core": "#f59e0b",
    "Catalytic site": "#f43f5e",
}


def _feature_bounds(feature: dict[str, Any]) -> tuple[int | None, int | None]:
    start = feature.get("start")
    end = feature.get("end")
    try:
        start_i = int(start) if start is not None else None
        end_i = int(end) if end is not None else start_i
    except (TypeError, ValueError):
        return None, None
    return start_i, end_i


def _extract_functional_annotations(
    feature_annotations: list[dict[str, Any]] | None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    if not feature_annotations:
        return [], []

    regions: list[dict[str, Any]] = []
    sites: list[dict[str, Any]] = []
    seen_regions: set[tuple[str, int, int]] = set()
    seen_sites: set[tuple[str, int, int]] = set()

    for feature in feature_annotations:
        feature_type = str(feature.get("type") or "").strip().lower()
        description = str(feature.get("description") or "").strip()
        search_text = f"{feature_type} {description}".lower()
        start, end = _feature_bounds(feature)
        if start is None:
            continue
        if end is None:
            end = start

        if feature_type == "active site" or "catalytic" in search_text:
            key = ("Catalytic site", start, end)
            if key in seen_sites:
                continue
            seen_sites.add(key)
            sites.append(
                {
                    "label": "Catalytic site",
                    "start": start,
                    "end": end,
                    "description": description or "Catalytic residue",
                    "color": _FUNCTIONAL_REGION_COLORS["Catalytic site"],
                }
            )
            continue

        if feature_type not in {"domain", "region"}:
            continue

        if "exonuclease" in search_text:
            label = "Exonuclease"
        elif "polymerase" in search_text:
            label = "Polymerase core"
        else:
            continue

        key = (label, start, end)
        if key in seen_regions:
            continue
        seen_regions.add(key)
        display_label = description or label
        regions.append(
            {
                "label": display_label,
                "start": start,
                "end": end,
                "description": description or label,
                "color": _FUNCTIONAL_REGION_COLORS[label],
            }
        )

    return regions, sites


def _coords_for_residue_range(
    pdb_coords: dict[int, dict[str, Any]],
    start: int,
    end: int,
) -> list[tuple[int, dict[str, Any]]]:
    lo, hi = sorted((start, end))
    return [
        (residue, pdb_coords[residue])
        for residue in sorted(pdb_coords)
        if lo <= residue <= hi
    ]


def _centroid(coords: list[tuple[int, dict[str, Any]]]) -> tuple[float, float, float]:
    if not coords:
        return 0.0, 0.0, 0.0
    count = float(len(coords))
    return (
        sum(entry["x"] for _, entry in coords) / count,
        sum(entry["y"] for _, entry in coords) / count,
        sum(entry["z"] for _, entry in coords) / count,
    )


def _add_functional_overlays(
    fig: go.Figure,
    pdb_coords: dict[int, dict[str, Any]],
    feature_annotations: list[dict[str, Any]] | None,
) -> None:
    regions, sites = _extract_functional_annotations(feature_annotations)

    for region in regions:
        coords = _coords_for_residue_range(
            pdb_coords,
            int(region["start"]),
            int(region["end"]),
        )
        if len(coords) < 2:
            continue

        fig.add_trace(
            go.Scatter3d(
                x=[coord["x"] for _, coord in coords],
                y=[coord["y"] for _, coord in coords],
                z=[coord["z"] for _, coord in coords],
                mode="lines",
                line={"color": region["color"], "width": 10},
                name=region["label"],
                showlegend=False,
                hovertemplate=(
                    f"<b>{region['label']}</b><br>"
                    f"Residues: {region['start']}-{region['end']}<br>"
                    f"{region['description']}"
                    "<extra></extra>"
                ),
            )
        )

        label_x, label_y, label_z = _centroid(coords)
        fig.add_trace(
            go.Scatter3d(
                x=[label_x],
                y=[label_y],
                z=[label_z],
                mode="text",
                text=[region["label"]],
                textposition="top center",
                textfont={"color": region["color"], "size": 12},
                name=f"{region['label']} label",
                showlegend=False,
                hoverinfo="skip",
            )
        )

    site_points: list[tuple[dict[str, Any], tuple[float, float, float]]] = []
    for site in sites:
        coords = _coords_for_residue_range(
            pdb_coords,
            int(site["start"]),
            int(site["end"]),
        )
        if not coords:
            continue
        site_points.append((site, _centroid(coords)))

    if not site_points:
        return

    fig.add_trace(
        go.Scatter3d(
            x=[point[1][0] for point in site_points],
            y=[point[1][1] for point in site_points],
            z=[point[1][2] for point in site_points],
            mode="markers",
            marker={
                "size": 9,
                "symbol": "cross",
                "color": _FUNCTIONAL_REGION_COLORS["Catalytic site"],
                "line": {"color": "#ffffff", "width": 1.5},
            },
            text=[
                (
                    f"{site['description']} "
                    f"({site['start']}"
                    + ("" if site["start"] == site["end"] else f"-{site['end']}")
                    + ")"
                )
                for site, _ in site_points
            ],
            name="Catalytic site",
            showlegend=False,
            hovertemplate="<b>%{text}</b><extra></extra>",
        )
    )

    label_x = sum(point[1][0] for point in site_points) / len(site_points)
    label_y = sum(point[1][1] for point in site_points) / len(site_points)
    label_z = sum(point[1][2] for point in site_points) / len(site_points)
    fig.add_trace(
        go.Scatter3d(
            x=[label_x],
            y=[label_y],
            z=[label_z],
            mode="text",
            text=["Catalytic site"],
            textposition="top center",
            textfont={"color": _FUNCTIONAL_REGION_COLORS["Catalytic site"], "size": 12},
            name="Catalytic site label",
            showlegend=False,
            hoverinfo="skip",
        )
    )


def build_3d_fingerprint(
    lineage_mutations: list[dict[str, Any]],
    backbone: list[dict[str, Any]],
    pdb_coords: dict[int, dict[str, Any]],
    variant_id: int,
    uniprot_id: str | None = None,
    structure_source: str = "Structure",
    feature_annotations: list[dict[str, Any]] | None = None,
) -> go.Figure:
    """
    Build a structure-first 3D mutation view using a C-alpha backbone trace.
    """
    fig = go.Figure()

    if backbone:
        fig.add_trace(
            go.Scatter3d(
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
                    "cmin": 0,
                    "cmax": 100,
                    "width": 5,
                },
                hoverinfo="skip",
                showlegend=False,
                name="C-alpha backbone",
            )
        )
        _add_functional_overlays(fig, pdb_coords, feature_annotations)

    generations = sorted({m["generation"] for m in lineage_mutations})
    mutation_type_styles = (
        ("non-synonymous", "Non-synonymous", 9, "diamond", "#f87171"),
        ("synonymous", "Synonymous", 6, "circle", "#94a3b8"),
    )
    for gen in generations:
        for mutation_type_key, mutation_type_label, marker_size, marker_symbol, marker_color in mutation_type_styles:
            xs: list[float] = []
            ys: list[float] = []
            zs: list[float] = []
            custom: list[list[Any]] = []

            for mutation in (m for m in lineage_mutations if m["generation"] == gen):
                trace_type = (
                    "synonymous"
                    if mutation["mutation_type"] == "synonymous"
                    else "non-synonymous"
                )
                if trace_type != mutation_type_key:
                    continue

                coord = pdb_coords.get(int(mutation["position"]))
                if coord is None:
                    continue

                xs.append(coord["x"])
                ys.append(coord["y"])
                zs.append(coord["z"])

                custom.append(
                    [
                        mutation.get("aa_change", f"{mutation['wt_aa']}{mutation['position']}{mutation['mut_aa']}"),
                        mutation["position"],
                        mutation["wt_aa"],
                        mutation["mut_aa"],
                        mutation["mutation_type"],
                        gen,
                        round(coord.get("plddt", 70.0), 1),
                    ]
                )

            if not xs:
                continue

            fig.add_trace(
                go.Scatter3d(
                    x=xs,
                    y=ys,
                    z=zs,
                    mode="markers",
                    name=f"Generation {gen} - {mutation_type_label}",
                    showlegend=False,
                    meta={
                        "trace_kind": "mutation",
                        "generation": gen,
                        "mutation_type": mutation_type_key,
                    },
                    marker={
                        "symbol": marker_symbol,
                        "size": marker_size,
                        "color": marker_color,
                        "opacity": 0.90,
                        "line": {"color": "rgba(0,0,0,0.30)", "width": 0.5},
                    },
                    customdata=custom,
                    hovertemplate=(
                        "<b>%{customdata[0]}</b><br>"
                        "Position: %{customdata[1]}<br>"
                        "Change: %{customdata[2]} -> %{customdata[3]}<br>"
                        "Type: %{customdata[4]}<br>"
                        "Generation: %{customdata[5]}<br>"
                        "pLDDT: %{customdata[6]}"
                        "<extra></extra>"
                    ),
                )
            )

    subtitle = f"UniProt {uniprot_id} · {structure_source}" if uniprot_id else structure_source
    fig.update_layout(
        title={
            "text": (
                f"Variant {variant_id} — 3D Mutation Fingerprint<br>"
                f"<sup>{subtitle} · "
                f"◆ non-synonymous &nbsp; ● synonymous &nbsp; "
                f"functional overlays for exonuclease / catalytic regions</sup>"
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
            "title": {"text": "Generation", "font": {"size": 11}},
            "bgcolor": "rgba(17,24,39,0.9)",
            "bordercolor": "#1e3050",
            "borderwidth": 1,
            "itemclick": "toggle",
            "itemdoubleclick": "toggleothers",
        },
        margin={"l": 0, "r": 0, "t": 80, "b": 0},
        height=780,
        hoverlabel={
            "bgcolor": "#162032",
            "bordercolor": "#1e3050",
            "font": {"family": "monospace", "size": 12, "color": "#d1dce8"},
        },
    )

    return fig


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    base_dir = Path(__file__).resolve().parent

    parser = argparse.ArgumentParser(
        description=(
            "Part 4a — Generate a circular mutation fingerprint HTML for a given variant. "
            "Sequence analysis is handled by mutation_analysis.py; "
            "this script owns only the Plotly visualisation."
        )
    )
    parser.add_argument(
        "--input",
        default=str(base_dir / "DE_BSU_Pol_Batch_1.tsv"),
        help="Variant table (.tsv or .json)",
    )
    parser.add_argument(
        "--variant-id", type=int, required=True,
        help="Variant ID to visualise",
    )
    parser.add_argument(
        "--wt-dna",
        default=str(base_dir / "pET-28a_BSU_DNA_Pol_I_WT.fa"),
        help="WT plasmid FASTA file",
    )
    parser.add_argument(
        "--wt-protein",
        default=str(base_dir / "O34996.fasta.txt"),
        help="WT protein FASTA file",
    )
    parser.add_argument(
        "--output", default=None,
        help="Output HTML file (default: variant_<id>_circular_fingerprint.html next to script)",
    )
    parser.add_argument(
        "--mutation-csv", default=None,
        help="Optional path to export a mutation CSV (produced by mutation_analysis.py)",
    )
    parser.add_argument(
        "--uniprot-id", default=None,
        help="UniProt accession for 3D structure fetch (auto-detected from WT FASTA header)",
    )
    parser.add_argument(
        "--pdb-file", default=None,
        help="Local PDB file for 3D mapping (skips AlphaFold / UniProt download)",
    )
    parser.add_argument(
        "--no-3d", action="store_true",
        help="Disable 3D mapping; render only the circular fingerprint",
    )
    parser.add_argument(
        "--linear", action="store_true",
        help=(
            "Render the linear (unfolded) fingerprint instead of the circular one. "
            "Mutations are colour-coded by generation across the variant lineage. "
            "Matches Figure 3 in the brief."
        ),
    )
    parser.add_argument(
        "--verbose", action="store_true",
        help="Enable verbose logging",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    input_path = Path(args.input).expanduser().resolve()
    wt_dna_path = Path(args.wt_dna).expanduser().resolve()
    wt_protein_path = Path(args.wt_protein).expanduser().resolve()

    # Default output filename reflects the chosen plot style
    plot_style = "linear_fingerprint" if args.linear else "circular_fingerprint"
    output_path = (
        Path(args.output).expanduser().resolve()
        if args.output
        else Path(__file__).resolve().parent
        / f"variant_{args.variant_id}_{plot_style}.html"
    )

    # ------------------------------------------------------------------
    # Step 1: Sequence analysis  (mutation_analysis.py)
    # ------------------------------------------------------------------
    _, wt_plasmid_seq = read_fasta(wt_dna_path)
    wt_header, wt_protein_seq = read_fasta(wt_protein_path)

    variants = load_variants_table(input_path)

    # ------------------------------------------------------------------
    # Step 2: Build figure
    # ------------------------------------------------------------------

    if args.linear:
        # Linear / unfolded view — colour mutations by generation (Figure 3)
        lineage_mutations = analyze_lineage_mutations(
            variants=variants,
            variant_id=int(args.variant_id),
            wt_plasmid_seq=wt_plasmid_seq,
            wt_protein_seq=wt_protein_seq,
        )
        # We still need the WT protein length for the backbone bar.
        analysis = analyze_target_variant(
            variants=variants,
            variant_id=int(args.variant_id),
            wt_plasmid_seq=wt_plasmid_seq,
            wt_protein_seq=wt_protein_seq,
        )
        fig = build_linear_fingerprint(
            lineage_mutations=lineage_mutations,
            wt_protein_len=len(analysis["wt_protein_sequence"]),
            variant_id=int(args.variant_id),
        )
        total_muts = len(lineage_mutations)
        nonsyn = sum(1 for m in lineage_mutations if m["mutation_type"] != "synonymous")
        syn = total_muts - nonsyn

    else:
        # Circular / folded view — original behaviour
        analysis = analyze_target_variant(
            variants=variants,
            variant_id=int(args.variant_id),
            wt_plasmid_seq=wt_plasmid_seq,
            wt_protein_seq=wt_protein_seq,
        )

        # Optionally export mutation CSV
        if args.mutation_csv:
            csv_path = Path(args.mutation_csv).expanduser().resolve()
            write_mutation_csv(csv_path, analysis)
            print(f"Mutation CSV: {csv_path}")

        uniprot_id = args.uniprot_id or infer_uniprot_id(wt_header)
        pdb_file = Path(args.pdb_file).expanduser().resolve() if args.pdb_file else None
        structure_coords, structure_info = resolve_structure(
            uniprot_id=uniprot_id,
            pdb_file=pdb_file,
            no_3d=bool(args.no_3d),
        )
        fig = build_plotly_figure(analysis, structure_coords, structure_info)
        total_muts = analysis["num_mutations"]
        nonsyn = analysis["num_nonsynonymous"]
        syn = analysis["num_synonymous"]

    # ------------------------------------------------------------------
    # Step 3: Save HTML
    # ------------------------------------------------------------------
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(output_path, include_plotlyjs="cdn", full_html=True)

    # Summary
    print(f"\nOutput:      {output_path}")
    print(f"Variant ID:  {args.variant_id}")
    print(f"Plot style:  {'linear (unfolded)' if args.linear else 'circular (folded)'}")
    print(
        f"Mutations:   total={total_muts}, "
        f"non-synonymous={nonsyn}, "
        f"synonymous={syn}"
    )
    if not args.linear:
        print(f"3D mapping:  {structure_info.get('status', 'n/a')}")
        if structure_info.get("source"):
            print(f"3D source:   {structure_info['source']}")


if __name__ == "__main__":
    main()