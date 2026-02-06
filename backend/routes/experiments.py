from flask import Blueprint, request, jsonify, session
from services.experiment_service import experiment_service

experiments_bp = Blueprint('experiments', __name__, url_prefix='/api/experiments')


def require_auth():
    """Check if user is authenticated"""
    user_id = session.get('user_id')
    if not user_id:
        return None
    return user_id


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
        
        # Build response with validation summary for frontend
        experiment_dict = experiment.to_dict(include_sequences=True)
        
        return jsonify({
            'success': True,
            'experiment': experiment_dict,
            'validation': {
                'isValid': experiment.validation_status == 'valid',
                'message': experiment.validation_message
            }
        }), 201
        
    except Exception as e:
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
        
    except Exception as e:
        return jsonify({'success': False, 'error': 'Server error'}), 500


@experiments_bp.route('/<experiment_id>', methods=['GET'])
def get_experiment(experiment_id: str):
    """Get a single experiment by ID"""
    user_id = require_auth()
    if not user_id:
        return jsonify({'success': False, 'error': 'Not authenticated'}), 401
    
    try:
        experiment = experiment_service.get_experiment_by_id(experiment_id, user_id)
        
        if not experiment:
            return jsonify({'success': False, 'error': 'Experiment not found'}), 404
        
        return jsonify({
            'success': True,
            'experiment': experiment.to_dict(include_sequences=True)
        }), 200
        
    except Exception as e:
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
        
    except Exception as e:
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
        
    except Exception as e:
        return jsonify({'success': False, 'error': 'Server error'}), 500
