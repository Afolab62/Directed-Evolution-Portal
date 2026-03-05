"""
Mutation fingerprint visualisation routes.

Routes registered here:
  GET  /api/experiments/<experiment_id>/fingerprint/<variant_id>
  GET  /api/experiments/<experiment_id>/fingerprint3d/<variant_id>
  GET  /api/experiments/<experiment_id>/fingerprint_linear/<variant_id>
"""
import json as _json
import uuid

from flask import request, jsonify, Response
from sqlalchemy.orm import joinedload

from database import db
from models.experiment import VariantData, safe_float
from services.experiment_service import experiment_service
from services.fingerprint_plot import resolve_structure, build_3d_fingerprint, build_linear_fingerprint
from services.uniprot_client import fetch_uniprot_features_json
from ._base import experiments_bp, require_auth


# ──────────────────────────────────────────────────────────────────────────────
# Shared helpers
# ──────────────────────────────────────────────────────────────────────────────

def _build_lineage(exp_uuid, var_uuid, extra_cols=None):
    """
    Walk the parent-chain to reconstruct the evolutionary lineage for a variant.

    Returns (selected_light, lineage_pvs) where lineage_pvs is a list of
    lightweight row-objects sorted by ascending generation.

    extra_cols: additional VariantData columns to include in the lightweight
    query (e.g. VariantData.protein_sequence).
    """
    base_cols = [
        VariantData.id,
        VariantData.plasmid_variant_index,
        VariantData.parent_plasmid_variant,
        VariantData.generation,
        VariantData.activity_score,
    ]
    if extra_cols:
        base_cols.extend(extra_cols)

    all_light = (
        db.query(*base_cols)
        .filter(VariantData.experiment_id == exp_uuid)
        .all()
    )

    by_pvi = {r.plasmid_variant_index: r for r in all_light}
    by_id  = {r.id: r for r in all_light}

    selected_light = by_id.get(var_uuid)
    if not selected_light:
        return None, None

    lineage_pvs: list = []
    current = selected_light
    visited: set = set()
    while current is not None and current.plasmid_variant_index not in visited:
        visited.add(current.plasmid_variant_index)
        lineage_pvs.append(current)
        parent_pvi = current.parent_plasmid_variant
        current = (
            by_pvi.get(parent_pvi)
            if parent_pvi is not None and parent_pvi >= 0
            else None
        )

    lineage_pvs.sort(key=lambda r: r.generation)
    return selected_light, lineage_pvs


def _load_lineage_with_mutations(lineage_pvs):
    """Load full ORM objects (with mutations eager-loaded) for the lineage."""
    lineage_ids = [r.id for r in lineage_pvs]
    lineage_variants = (
        db.query(VariantData)
        .filter(VariantData.id.in_(lineage_ids))
        .options(joinedload(VariantData.mutations))
        .all()
    )
    lineage_by_id = {v.id: v for v in lineage_variants}
    return [lineage_by_id[r.id] for r in lineage_pvs if r.id in lineage_by_id]


def _delta_walk_nonsynonymous(lineage):
    """
    Delta walk that tracks only non-synonymous mutations.
    Returns a list of fingerprint dicts with generationIntroduced.
    """
    fingerprint: list = []
    seen_keys: set = set()
    prev_mut_keys: set = set()

    for step in lineage:
        step_mut_keys = {
            (m.position, m.wild_type, m.mutant)
            for m in (step.mutations or [])
            if m.mutation_type != 'synonymous'
        }
        new_keys = step_mut_keys - prev_mut_keys

        for m in (step.mutations or []):
            if m.mutation_type == 'synonymous':
                continue
            key = (m.position, m.wild_type, m.mutant)
            if key in new_keys and key not in seen_keys:
                seen_keys.add(key)
                fingerprint.append({
                    'position': m.position,
                    'wildType': m.wild_type,
                    'mutant': m.mutant,
                    'mutationType': m.mutation_type,
                    'generationIntroduced': step.generation,
                })

        prev_mut_keys = step_mut_keys

    fingerprint.sort(key=lambda m: m['position'])
    return fingerprint


def _delta_walk_all(lineage):
    """
    Delta walk that collects ALL mutations (synonymous + non-synonymous).
    Returns a list of mutation dicts including codon info and generation.
    """
    all_mutations_raw: dict = {}
    prev_keys: set = set()

    for step in lineage:
        step_keys = {
            (m.position, m.wild_type, m.mutant, m.mutation_type)
            for m in (step.mutations or [])
        }
        new_keys = step_keys - prev_keys
        for m in (step.mutations or []):
            key = (m.position, m.wild_type, m.mutant, m.mutation_type)
            if key in new_keys and key not in all_mutations_raw:
                all_mutations_raw[key] = {
                    'position':      m.position,
                    'wt_aa':         m.wild_type,
                    'mut_aa':        m.mutant,
                    'mutation_type': m.mutation_type,
                    'wt_codon':      m.wt_codon or '',
                    'mut_codon':     m.mut_codon or '',
                    'aa_change':     f"{m.wild_type}{m.position}{m.mutant}",
                    'generation':    step.generation,
                }
        prev_keys = step_keys

    return sorted(all_mutations_raw.values(), key=lambda m: m['position'])


# ──────────────────────────────────────────────────────────────────────────────
# Routes
# ──────────────────────────────────────────────────────────────────────────────

@experiments_bp.route('/<experiment_id>/fingerprint/<variant_id>', methods=['GET'])
def get_mutation_fingerprint(experiment_id: str, variant_id: str):
    """
    Return the mutation fingerprint for a specific variant.
    Reconstructs the evolutionary lineage and assigns generation_introduced
    for each non-synonymous mutation.
    """
    user_id = require_auth()
    if not user_id:
        return jsonify({'success': False, 'error': 'Not authenticated'}), 401

    try:
        experiment = experiment_service.get_experiment_by_id(experiment_id, user_id)
        if not experiment:
            return jsonify({'success': False, 'error': 'Experiment not found'}), 404

        exp_uuid = uuid.UUID(experiment_id)
        var_uuid = uuid.UUID(variant_id)

        selected_light, lineage_pvs = _build_lineage(exp_uuid, var_uuid)
        if selected_light is None:
            return jsonify({'success': False, 'error': 'Variant not found'}), 404

        lineage = _load_lineage_with_mutations(lineage_pvs)
        fingerprint = _delta_walk_nonsynonymous(lineage)

        protein_length = (
            len(experiment.wt_protein_sequence.strip())
            if experiment.wt_protein_sequence else 0
        )

        return jsonify({
            'success': True,
            'variantId': str(selected_light.id),
            'plasmidVariantIndex': selected_light.plasmid_variant_index,
            'generation': selected_light.generation,
            'activityScore': safe_float(selected_light.activity_score),
            'proteinLength': protein_length,
            'lineageLength': len(lineage),
            'fingerprint': fingerprint,
        }), 200

    except Exception:
        import traceback
        traceback.print_exc()
        return jsonify({'success': False, 'error': 'Server error'}), 500


@experiments_bp.route('/<experiment_id>/fingerprint3d/<variant_id>', methods=['GET'])
def get_mutation_fingerprint_3d(experiment_id: str, variant_id: str):
    """
    Return a Plotly JSON figure for the circular mutation fingerprint with an
    optional 3D structure panel (AlphaFold / PDB via protein_accession).
    """
    user_id = require_auth()
    if not user_id:
        return jsonify({'success': False, 'error': 'Not authenticated'}), 401

    try:
        experiment = experiment_service.get_experiment_by_id(experiment_id, user_id)
        if not experiment:
            return jsonify({'success': False, 'error': 'Experiment not found'}), 404

        exp_uuid = uuid.UUID(experiment_id)
        var_uuid = uuid.UUID(variant_id)

        selected_light, lineage_pvs = _build_lineage(
            exp_uuid, var_uuid,
            extra_cols=[VariantData.protein_sequence]
        )
        if selected_light is None:
            return jsonify({'success': False, 'error': 'Variant not found'}), 404

        lineage = _load_lineage_with_mutations(lineage_pvs)
        mutations = _delta_walk_all(lineage)

        num_synonymous    = sum(1 for m in mutations if m['mutation_type'] == 'synonymous')
        num_nonsynonymous = len(mutations) - num_synonymous

        # Resolve 3-D structure
        uniprot_id = experiment.protein_accession or None
        no_3d = request.args.get('no_3d', 'false').lower() == 'true'
        try:
            backbone, pdb_coords, structure_info = resolve_structure(
                uniprot_id=uniprot_id,
                no_3d=no_3d,
            )
        except Exception as ex:
            import logging
            logging.getLogger(__name__).warning("Structure resolution failed: %s", ex)
            backbone, pdb_coords, structure_info = [], {}, {'status': f'3D mapping failed: {ex}'}

        # UniProt feature annotations
        feature_annotations: list = []
        if uniprot_id:
            try:
                feature_annotations = fetch_uniprot_features_json(uniprot_id)
            except Exception:
                pass

        # Build Plotly figure
        variant_pvi = selected_light.plasmid_variant_index
        highlight_position = request.args.get('highlight', type=int)
        fig = build_3d_fingerprint(
            lineage_mutations=mutations,
            backbone=backbone,
            pdb_coords=pdb_coords,
            variant_id=variant_pvi,
            uniprot_id=uniprot_id,
            structure_source=structure_info.get('source', structure_info.get('status', 'Structure')),
            feature_annotations=feature_annotations,
            highlight_position=highlight_position,
        )

        if request.args.get('format') == 'html':
            html = fig.to_html(include_plotlyjs='cdn', full_html=True, config={
                'responsive': True,
                'modeBarButtonsToRemove': ['select2d', 'lasso2d'],
            })
            return Response(html, mimetype='text/html')

        fig_json = _json.loads(fig.to_json())

        return jsonify({
            'success': True,
            'figure': fig_json,
            'variantId': str(selected_light.id),
            'plasmidVariantIndex': selected_light.plasmid_variant_index,
            'generation': selected_light.generation,
            'activityScore': safe_float(selected_light.activity_score),
            'numMutations': len(mutations),
            'numSynonymous': num_synonymous,
            'numNonSynonymous': num_nonsynonymous,
            'structureStatus': structure_info.get('status', ''),
            'lineageLength': len(lineage),
        }), 200

    except Exception:
        import traceback
        traceback.print_exc()
        return jsonify({'success': False, 'error': 'Server error'}), 500


@experiments_bp.route('/<experiment_id>/fingerprint_linear/<variant_id>', methods=['GET'])
def get_mutation_fingerprint_linear(experiment_id: str, variant_id: str):
    """
    Return a Plotly JSON figure for the linear mutation fingerprint.
    Triangles on a horizontal backbone, one colour per generation.
    Supports ?format=html for a self-contained interactive HTML page.
    Clicking a mutation triangle posts a message to the parent frame so React
    can auto-switch to the 3D tab with the selected residue highlighted.
    """
    user_id = require_auth()
    if not user_id:
        return jsonify({'success': False, 'error': 'Not authenticated'}), 401

    try:
        experiment = experiment_service.get_experiment_by_id(experiment_id, user_id)
        if not experiment:
            return jsonify({'success': False, 'error': 'Experiment not found'}), 404

        exp_uuid = uuid.UUID(experiment_id)
        var_uuid = uuid.UUID(variant_id)

        wt_protein_len = len((experiment.wt_protein_sequence or '').strip()) or 500

        selected_light, lineage_pvs = _build_lineage(exp_uuid, var_uuid)
        if selected_light is None:
            return jsonify({'success': False, 'error': 'Variant not found'}), 404

        lineage = _load_lineage_with_mutations(lineage_pvs)
        mutations = _delta_walk_all(lineage)

        _win_start = request.args.get('window_start', type=int)
        _win_end   = request.args.get('window_end',   type=int)

        fig = build_linear_fingerprint(
            lineage_mutations=mutations,
            wt_protein_len=wt_protein_len,
            variant_id=selected_light.plasmid_variant_index,
            window_start=_win_start,
            window_end=_win_end,
        )

        if request.args.get('format') == 'html':
            html = fig.to_html(include_plotlyjs='cdn', full_html=True, config={
                'responsive': True,
                'modeBarButtonsToRemove': ['select2d', 'lasso2d'],
            })
            # Inject a click handler that fires window.parent.postMessage so
            # the parent React component can auto-switch to the 3D tab.
            _click_script = """<script>
(function () {
  function attach() {
    var divs = document.querySelectorAll('.js-plotly-plot');
    if (!divs.length) { setTimeout(attach, 150); return; }
    divs.forEach(function (div) {
      div.on('plotly_click', function (ev) {
        if (!ev || !ev.points || !ev.points.length) return;
        var pt = ev.points[0];
        var pos = pt.x != null ? Math.round(pt.x) : null;
        var label = pt.hovertext || null;
        if (pos != null) {
          window.parent.postMessage(
            { type: 'dem_mutation_click', position: pos, label: label },
            '*'
          );
        }
      });
    });
  }
  if (document.readyState === 'complete') { attach(); }
  else { window.addEventListener('load', attach); }
})();
</script>"""
            html = html.replace('</body>', _click_script + '\n</body>')
            return Response(html, mimetype='text/html')

        fig_json = _json.loads(fig.to_json())

        return jsonify({
            'success': True,
            'figure': fig_json,
            'plasmidVariantIndex': selected_light.plasmid_variant_index,
            'generation': selected_light.generation,
            'activityScore': safe_float(selected_light.activity_score),
            'numMutations': len(mutations),
            'wtProteinLen': wt_protein_len,
        }), 200

    except Exception:
        import traceback
        traceback.print_exc()
        return jsonify({'success': False, 'error': 'Server error'}), 500
