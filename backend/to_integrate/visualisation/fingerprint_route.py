"""
app/routes/fingerprint.py
--------------------------
Part 4a — Fingerprint Blueprint
=================================

Two routes:

  GET /fingerprint/
      Lists all variants with generation + activity score.
      User clicks a row → opens the 3D plot for that variant.

  GET /fingerprint/<int:variant_id>
      Runs the full lineage analysis for one variant, fetches
      the AlphaFold PDB for the experiment's UniProt ID (auto-
      downloaded and cached), and renders the 3D plot page.

Where does this live in the app?
---------------------------------
The app already has a fingerprint_bp registered in __init__.py:

    from .routes.fingerprint import fingerprint_bp
    app.register_blueprint(fingerprint_bp)

This file IS that blueprint.  It slots directly after the staging /
main pipeline — once variants are in the DB the user navigates to
/fingerprint/ to pick one and see its 3D mutation map.

Does the plot auto-fetch from AlphaFold?
-----------------------------------------
Yes. resolve_structure() in fingerprint_plot.py tries:
  1. AlphaFold EBI API  (AF-{UNIPROT_ID}-F1-model_v4.pdb, then v3/v2/v1)
  2. UniProt PDB cross-reference as fallback
The downloaded PDB is cached in STRUCTURE_CACHE_DIR (see config.py)
so subsequent requests for the same protein are instant.

The UniProt ID comes from the staging record stored in the DB at
the end of the plasmid-validation step.
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

from flask import Blueprint, render_template, abort, current_app, jsonify, request

# ── Analysis layer (same package) ──────────────────────────────────────────
from app.analysis.mutation_analysis import (
    analyze_target_variant,
    assign_generation_to_target_mutations,
    get_variant_lineage,
)
from app.analysis.fingerprint_plot import (
    resolve_structure,
    build_3d_fingerprint,
)

# ── DB helpers ──────────────────────────────────────────────────────────────
from app.services.fingerprint_db import (
    load_variants_from_db,
    get_all_variants_with_scores,
    get_staging_info,
)
from app.services.staging import fetch_uniprot_features_json

logger = logging.getLogger(__name__)

fingerprint_bp = Blueprint("fingerprint", __name__, url_prefix="/fingerprint")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _parse_pdb_for_plot(pdb_path: Path) -> tuple[list[dict], dict[int, dict]]:
    """
    Parse a PDB file and return:
      backbone  — list of {r, x, y, z, plddt} dicts for the Cα ribbon
      coords    — {residue_number: {x, y, z, plddt}} for mutation mapping
    """
    coords: dict[int, dict] = {}
    with open(pdb_path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            try:
                res_num = int(line[22:26].strip())
                coords[res_num] = {
                    "x":     round(float(line[30:38]), 3),
                    "y":     round(float(line[38:46]), 3),
                    "z":     round(float(line[46:54]), 3),
                    "plddt": round(float(line[60:66]), 1),
                }
            except (ValueError, IndexError):
                continue

    backbone = [
        {"r": r, **c}
        for r, c in sorted(coords.items())
    ]
    return backbone, coords


def _attach_3d_coords(mutations: list[dict], coords: dict[int, dict]) -> list[dict]:
    """Add x/y/z/plddt keys to each mutation dict from the PDB coords."""
    for m in mutations:
        c = coords.get(m["position"])
        m["x"]     = c["x"]     if c else None
        m["y"]     = c["y"]     if c else None
        m["z"]     = c["z"]     if c else None
        m["plddt"] = c["plddt"] if c else None
    return mutations


# ---------------------------------------------------------------------------
# GET /fingerprint/  — variant selector
# ---------------------------------------------------------------------------

@fingerprint_bp.route("/")
def index():
    """
    Show every variant with its generation, parent ID and activity score.
    The user clicks a row to open the 3D fingerprint for that variant.
    """
    db_path = Path(current_app.config["DB_PATH"])
    variants = get_all_variants_with_scores(db_path)
    requested_uniprot = (request.args.get("uniprot_id") or "").strip().upper()
    staging = get_staging_info(db_path)
    return render_template(
        "fingerprint/index.html",
        variants=variants,
        uniprot_id=requested_uniprot or (staging or {}).get("uniprot_id", ""),
    )


# ---------------------------------------------------------------------------
# GET /fingerprint/<variant_id>  — 3D mutation plot
# ---------------------------------------------------------------------------

@fingerprint_bp.route("/<int:variant_id>")
def show(variant_id: int):
    """
    Full pipeline for one variant:

    1. Load all variants from DB (needed for lineage walk).
    2. Load WT sequences + UniProt ID from the staging record.
    3. Auto-fetch AlphaFold PDB for the UniProt ID (cached on disk).
    4. Parse Cα backbone + per-residue pLDDT from the PDB.
    5. Run analyze_lineage_mutations() to get mutations tagged with the
       generation they were *first introduced* (not just present in 295).
    6. Attach 3D coordinates to each mutation.
    7. Render template with pre-serialised JSON — no extra AJAX needed.
    """
    db_path = Path(current_app.config["DB_PATH"])

    # ── 1. Load variants ──────────────────────────────────────────────────
    variants = load_variants_from_db(db_path)
    if variant_id not in variants:
        abort(404, description=f"Variant {variant_id} not found.")

    # ── 2. WT sequences + UniProt ID from staging record ─────────────────
    staging = get_staging_info(db_path)
    if not staging:
        abort(500, description="No staging record found — run the staging step first.")

    wt_protein_seq = staging["wt_protein_seq"]
    wt_plasmid_seq = staging["wt_plasmid_seq"]
    requested_uniprot = (request.args.get("uniprot_id") or "").strip().upper()
    uniprot_id = requested_uniprot or staging["uniprot_id"]
    feature_annotations: list[dict] = []
    if uniprot_id:
        try:
            feature_annotations = fetch_uniprot_features_json(uniprot_id)
        except Exception as exc:
            logger.warning("Could not load UniProt features for %s: %s", uniprot_id, exc)

    # ── 3. Structure lookup — auto-fetched and cached ────────────────────
    # resolve_structure() uses the configured structure source first, then
    # falls back to UniProt-linked PDB entries when available.
    # Returns (coords_dict, info_dict) where coords_dict is None on failure.
    structure_coords_raw, structure_info = resolve_structure(
        uniprot_id=uniprot_id,
        pdb_file=None,
        no_3d=False,
    )

    # ── 4. Parse backbone + per-residue data from the cached PDB ─────────
    backbone: list[dict] = []
    pdb_coords: dict[int, dict] = {}

    if structure_info.get("pdb_path"):
        # fingerprint_plot.resolve_structure can return the path in info
        backbone, pdb_coords = _parse_pdb_for_plot(
            Path(structure_info["pdb_path"])
        )
    elif structure_coords_raw:
        # Fallback: coords_raw is {res: (x,y,z)} — no pLDDT available
        for res, (x, y, z) in sorted(structure_coords_raw.items()):
            entry = {"x": round(x, 3), "y": round(y, 3), "z": round(z, 3), "plddt": 70.0}
            pdb_coords[res] = entry
            backbone.append({"r": res, **entry})

    # ── 5. Final target mutations + first-seen generation ─────────────────
    # Use codon-level target mutations so synonymous changes survive into the
    # 3D plot, then annotate each final mutation with the earliest lineage
    # generation where that final codon state appears.
    try:
        target_analysis = analyze_target_variant(
            variants=variants,
            variant_id=variant_id,
            wt_plasmid_seq=wt_plasmid_seq,
            wt_protein_seq=wt_protein_seq,
        )
        mutations = assign_generation_to_target_mutations(
            variants=variants,
            variant_id=variant_id,
            wt_plasmid_seq=wt_plasmid_seq,
            wt_protein_seq=wt_protein_seq,
            target_mutations=target_analysis["mutations"],
            default_generation=target_analysis["generation"],
        )
    except Exception as exc:
        logger.error("Target mutation analysis failed for variant %d: %s", variant_id, exc)
        abort(500, description=f"Target mutation analysis failed: {exc}")

    # ── 6. Attach 3D coordinates ──────────────────────────────────────────
    mutations = _attach_3d_coords(mutations, pdb_coords)

    # ── 7. Build the Plotly figure using build_3d_fingerprint() ──────────
    fig = build_3d_fingerprint(
        lineage_mutations=mutations,
        backbone=backbone,
        pdb_coords=pdb_coords,
        variant_id=variant_id,
        uniprot_id=uniprot_id,
        structure_source=structure_info.get("source", "Structure"),
        feature_annotations=feature_annotations,
    )
    plot_json = fig.to_json()   # passed into the template as a safe JSON string

    # ── 8. Pre-compute generation counts for the sidebar ─────────────────
    gen_counts = {g: 0 for g in range(1, 11)}
    for m in mutations:
        g = m.get("generation", 0)
        if 1 <= g <= 10:
            gen_counts[g] += 1

    var_info = variants[variant_id]
    lineage  = get_variant_lineage(variants, variant_id)

    return render_template(
        "fingerprint/show.html",
        # Template context
        variant_id       = variant_id,
        generation       = var_info["generation"],
        parent_id        = var_info["parent_variant_id"],
        uniprot_id       = uniprot_id,
        structure_source = structure_info.get("source", "Unknown"),
        lineage          = lineage,
        num_mutations    = target_analysis["num_mutations"],
        num_nonsynonymous= target_analysis["num_nonsynonymous"],
        num_synonymous   = target_analysis["num_synonymous"],
        # Plotly figure JSON — rendered by Plotly.newPlot in the template
        plot_json        = plot_json,
        # Raw data also passed for the sidebar generation legend
        gen_counts_json  = json.dumps(gen_counts),
    )


# ---------------------------------------------------------------------------
# GET /fingerprint/<variant_id>/data  — raw JSON (for API / future AJAX)
# ---------------------------------------------------------------------------

@fingerprint_bp.route("/<int:variant_id>/data")
def data(variant_id: int):
    """Return the full mutation + backbone payload as JSON."""
    db_path   = Path(current_app.config["DB_PATH"])

    variants = load_variants_from_db(db_path)
    if variant_id not in variants:
        return jsonify(error=f"Variant {variant_id} not found"), 404

    staging = get_staging_info(db_path)
    if not staging:
        return jsonify(error="No staging record"), 500

    structure_coords_raw, structure_info = resolve_structure(
        uniprot_id=staging["uniprot_id"],
        pdb_file=None,
        no_3d=False,
    )

    backbone, pdb_coords = [], {}
    if structure_info.get("pdb_path"):
        backbone, pdb_coords = _parse_pdb_for_plot(Path(structure_info["pdb_path"]))
    elif structure_coords_raw:
        for res, (x, y, z) in sorted(structure_coords_raw.items()):
            entry = {"x": round(x,3), "y": round(y,3), "z": round(z,3), "plddt": 70.0}
            pdb_coords[res] = entry
            backbone.append({"r": res, **entry})

    target = analyze_target_variant(
        variants=variants,
        variant_id=variant_id,
        wt_plasmid_seq=staging["wt_plasmid_seq"],
        wt_protein_seq=staging["wt_protein_seq"],
    )
    mutations = assign_generation_to_target_mutations(
        variants=variants,
        variant_id=variant_id,
        wt_plasmid_seq=staging["wt_plasmid_seq"],
        wt_protein_seq=staging["wt_protein_seq"],
        target_mutations=target["mutations"],
        default_generation=target["generation"],
    )
    mutations = _attach_3d_coords(mutations, pdb_coords)

    return jsonify(
        variant_id       = variant_id,
        generation       = variants[variant_id]["generation"],
        uniprot_id       = staging["uniprot_id"],
        structure        = structure_info,
        num_mutations    = target["num_mutations"],
        num_nonsynonymous= target["num_nonsynonymous"],
        num_synonymous   = target["num_synonymous"],
        mutations        = mutations,
        backbone         = backbone,
    )