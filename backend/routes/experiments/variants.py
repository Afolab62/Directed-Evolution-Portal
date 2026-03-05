"""
Variant access routes.

Routes registered here:
  GET  /api/experiments/<experiment_id>/variants
  GET  /api/experiments/<experiment_id>/top-performers
"""
from flask import request, jsonify
import uuid

from database import db
from models.experiment import VariantData
from services.experiment_service import experiment_service
from sqlalchemy.orm import joinedload
from ._base import experiments_bp, require_auth


@experiments_bp.route('/<experiment_id>/variants', methods=['GET'])
def get_experiment_variants(experiment_id: str):
    """Get all variant data for an experiment (paginated)"""
    user_id = require_auth()
    if not user_id:
        return jsonify({'success': False, 'error': 'Not authenticated'}), 401

    try:
        experiment = experiment_service.get_experiment_by_id(experiment_id, user_id)
        if not experiment:
            return jsonify({'success': False, 'error': 'Experiment not found'}), 404

        exp_uuid = uuid.UUID(experiment_id) if isinstance(experiment_id, str) else experiment_id
        limit = min(int(request.args.get('limit', 1000)), 5000)
        include_mutations = request.args.get('include_mutations', 'false').lower() == 'true'

        variants = db.query(VariantData).filter_by(
            experiment_id=exp_uuid
        ).order_by(
            VariantData.generation.asc(),
            VariantData.activity_score.desc().nullslast()
        ).limit(limit).all()

        return jsonify({
            'success': True,
            'variants': [v.to_dict(include_mutations=include_mutations) for v in variants]
        }), 200

    except Exception:
        return jsonify({'success': False, 'error': 'Server error'}), 500


@experiments_bp.route('/<experiment_id>/top-performers', methods=['GET'])
def get_top_performers(experiment_id: str):
    """Get top performing variants by activity score"""
    user_id = require_auth()
    if not user_id:
        return jsonify({'success': False, 'error': 'Not authenticated'}), 401

    try:
        experiment = experiment_service.get_experiment_by_id(experiment_id, user_id)
        if not experiment:
            return jsonify({'success': False, 'error': 'Experiment not found'}), 404

        limit = min(int(request.args.get('limit', 10)), 50)
        include_mutations = request.args.get('include_mutations', 'true').lower() == 'true'

        exp_uuid = uuid.UUID(experiment_id) if isinstance(experiment_id, str) else experiment_id

        query = db.query(VariantData)
        if include_mutations:
            query = query.options(joinedload(VariantData.mutations))

        variants = query.filter_by(
            experiment_id=exp_uuid,
            is_control=False
        ).filter(
            VariantData.activity_score.isnot(None)
        ).order_by(
            VariantData.activity_score.desc()
        ).limit(limit).all()

        for v in variants:
            mut_list = v.mutations if include_mutations else []
            print(f"[TOP] variant pvi={v.plasmid_variant_index} gen={v.generation} "
                  f"mutations={len(mut_list)} activity={v.activity_score:.3f}")

        return jsonify({
            'success': True,
            'topPerformers': [
                v.to_dict(include_sequences=True, include_mutations=include_mutations)
                for v in variants
            ]
        }), 200

    except Exception:
        return jsonify({'success': False, 'error': 'Server error'}), 500
