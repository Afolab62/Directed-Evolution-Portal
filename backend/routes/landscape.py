import json as _json
import uuid
from flask import Blueprint, request, jsonify, session
from database import db
from models.experiment import VariantData
from services.landscape_service import build_landscape_figure

landscape_bp = Blueprint('landscape', __name__, url_prefix='/api/experiments')


def require_auth():
    user_id = session.get('user_id')
    if not user_id:
        return None
    return user_id


@landscape_bp.route('/<experiment_id>/landscape', methods=['GET'])
def get_fitness_landscape(experiment_id: str):
    """
    Compute and return a complete Plotly figure JSON for the 3D activity
    landscape.  The figure includes per-generation animation frames, three
    z-mode transforms (robust / raw / normalized), play buttons, and a
    generation slider — matching the output of 3d_landscape.py.

    Query params:
        method: "pca" (default) | "tsne" | "umap"
    """
    user_id = require_auth()
    if not user_id:
        return jsonify({'success': False, 'error': 'Not authenticated'}), 401

    try:
        method = request.args.get('method', 'pca').lower()
        if method not in ('pca', 'tsne', 'umap'):
            return jsonify({'success': False, 'error': 'method must be pca, tsne, or umap'}), 400

        exp_uuid = uuid.UUID(experiment_id)

        variants = (
            db.query(VariantData)
            .filter_by(experiment_id=exp_uuid, qc_status='passed', is_control=False)
            .filter(VariantData.protein_sequence.isnot(None))
            .filter(VariantData.activity_score.isnot(None))
            .all()
        )

        if len(variants) < 3:
            return jsonify({
                'success': False,
                'error': f'Need at least 3 analysed variants (found {len(variants)}). '
                         'Run sequence analysis first.'
            }), 422

        fig = build_landscape_figure(
            sequences      =[v.protein_sequence for v in variants],
            activity_scores=[float(v.activity_score) for v in variants],
            generations    =[v.generation for v in variants],
            variant_indices=[v.plasmid_variant_index for v in variants],
            method=method,
        )

        fig_json = _json.loads(fig.to_json())
        return jsonify({
            'success':       True,
            'figure':        fig_json,
            'variant_count': len(variants),
            'method':        method,
        }), 200

    except ValueError as e:
        return jsonify({'success': False, 'error': str(e)}), 422
    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({'success': False, 'error': f'Landscape computation failed: {str(e)}'}), 500
