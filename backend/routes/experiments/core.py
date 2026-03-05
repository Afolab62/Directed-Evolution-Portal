"""
Core CRUD routes for experiments.

Routes registered here:
  POST   /api/experiments
  GET    /api/experiments
  GET    /api/experiments/<experiment_id>
  PATCH  /api/experiments/<experiment_id>
  DELETE /api/experiments/<experiment_id>
"""
from flask import request, jsonify
import uuid

from database import db
from models.experiment import VariantData
from services.experiment_service import experiment_service
from ._base import experiments_bp, require_auth


@experiments_bp.route('', methods=['POST'])
def create_experiment():
    """Create a new experiment"""
    user_id = require_auth()
    if not user_id:
        return jsonify({'success': False, 'error': 'Not authenticated'}), 401

    try:
        data = request.get_json()

        name = data.get('name', '').strip()
        protein_accession = data.get('proteinAccession', '').strip()
        plasmid_sequence = data.get('plasmidSequence', '').strip()
        plasmid_name = data.get('plasmidName', '').strip() or None
        fetch_features = data.get('fetchFeatures', True)

        if not name or not protein_accession or not plasmid_sequence:
            return jsonify({
                'success': False,
                'error': 'Name, protein accession, and plasmid sequence are required'
            }), 400

        # ── Validate FASTA structure before hitting UniProt ──────────────────
        first_non_empty = next(
            (ln.strip() for ln in plasmid_sequence.splitlines() if ln.strip()), ''
        )
        if not first_non_empty.startswith('>'):
            return jsonify({
                'success': False,
                'error': (
                    'Plasmid file must be in FASTA format '
                    '(first line must start with ">"). '
                    'If you have a raw DNA sequence, add a header line such as ">MyPlasmid" before it.'
                )
            }), 400

        dna_chars = set()
        for ln in plasmid_sequence.splitlines():
            if not ln.strip().startswith('>'):
                dna_chars.update(ln.strip().upper())
        invalid_chars = dna_chars - set('ACGTN\r\n ')
        if invalid_chars:
            return jsonify({
                'success': False,
                'error': (
                    f'Plasmid sequence contains invalid DNA characters: '
                    f'{", ".join(sorted(invalid_chars))}. '
                    'DNA sequences may only contain A, C, G, T and N.'
                )
            }), 400

        experiment, error = experiment_service.create_experiment(
            user_id=user_id,
            name=name,
            protein_accession=protein_accession,
            plasmid_sequence=plasmid_sequence,
            plasmid_name=plasmid_name,
            fetch_features=fetch_features
        )

        if error:
            return jsonify({'success': False, 'error': error}), 400

        experiment_dict = experiment.to_dict(include_sequences=True)

        return jsonify({
            'success': True,
            'experiment': experiment_dict,
            'validation': {
                'isValid': experiment.validation_status == 'valid',
                'message': experiment.validation_message
            }
        }), 201

    except Exception:
        return jsonify({'success': False, 'error': 'Server error'}), 500


@experiments_bp.route('', methods=['GET'])
def list_experiments():
    """Get all experiments for the current user"""
    user_id = require_auth()
    if not user_id:
        return jsonify({'success': False, 'error': 'Not authenticated'}), 401

    try:
        limit = min(int(request.args.get('limit', 100)), 100)
        offset = int(request.args.get('offset', 0))

        experiments = experiment_service.get_experiments_by_user(user_id, limit=limit, offset=offset)

        return jsonify({
            'success': True,
            'experiments': [exp.to_dict() for exp in experiments]
        }), 200

    except Exception:
        return jsonify({'success': False, 'error': 'Server error'}), 500


@experiments_bp.route('/<experiment_id>', methods=['GET'])
def get_experiment(experiment_id: str):
    """Get a single experiment by ID with its variants"""
    user_id = require_auth()
    if not user_id:
        return jsonify({'success': False, 'error': 'Not authenticated'}), 401

    try:
        experiment = experiment_service.get_experiment_by_id(experiment_id, user_id)

        if not experiment:
            return jsonify({'success': False, 'error': 'Experiment not found'}), 404

        include_variants = request.args.get('include_variants', 'true').lower() != 'false'

        variants = []
        if include_variants:
            try:
                exp_uuid = uuid.UUID(experiment_id) if isinstance(experiment_id, str) else experiment_id
                limit = min(int(request.args.get('limit', 1000)), 5000)

                variants = db.query(VariantData).filter_by(
                    experiment_id=exp_uuid
                ).order_by(
                    VariantData.generation.asc(),
                    VariantData.activity_score.desc().nullslast()
                ).limit(limit).all()

                print(f"Loaded {len(variants)} variants for experiment {experiment_id}")
            except Exception as e:
                print(f"Error fetching variants: {e}")
                variants = []

        return jsonify({
            'success': True,
            'experiment': experiment.to_dict(include_sequences=True),
            'variants': [v.to_dict(include_mutations=False) for v in variants]
        }), 200

    except Exception as e:
        print(f"Error in get_experiment: {e}")
        import traceback
        traceback.print_exc()
        return jsonify({'success': False, 'error': 'Server error'}), 500


@experiments_bp.route('/<experiment_id>', methods=['PATCH'])
def update_experiment(experiment_id: str):
    """Update experiment metadata"""
    user_id = require_auth()
    if not user_id:
        return jsonify({'success': False, 'error': 'Not authenticated'}), 401

    try:
        data = request.get_json()

        name = data.get('name')
        plasmid_name = data.get('plasmidName')

        experiment = experiment_service.update_experiment(
            experiment_id=experiment_id,
            user_id=user_id,
            name=name,
            plasmid_name=plasmid_name
        )

        if not experiment:
            return jsonify({'success': False, 'error': 'Experiment not found'}), 404

        return jsonify({
            'success': True,
            'experiment': experiment.to_dict(include_sequences=True)
        }), 200

    except Exception:
        return jsonify({'success': False, 'error': 'Server error'}), 500


@experiments_bp.route('/<experiment_id>', methods=['DELETE'])
def delete_experiment(experiment_id: str):
    """Delete an experiment"""
    user_id = require_auth()
    if not user_id:
        return jsonify({'success': False, 'error': 'Not authenticated'}), 401

    try:
        success = experiment_service.delete_experiment(experiment_id, user_id)

        if not success:
            return jsonify({'success': False, 'error': 'Experiment not found'}), 404

        return jsonify({'success': True, 'message': 'Experiment deleted'}), 200

    except Exception:
        return jsonify({'success': False, 'error': 'Server error'}), 500
