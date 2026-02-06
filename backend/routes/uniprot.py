from flask import Blueprint, jsonify
from services.uniprot_client import (
    fetch_uniprot_protein_metadata,
    UniProtNotFound,
    UniProtNetworkError
)

uniprot_bp = Blueprint('uniprot', __name__, url_prefix='/api/uniprot')


@uniprot_bp.route('/<accession>', methods=['GET'])
def get_protein(accession: str):
    """Fetch protein data from UniProt"""
    if not accession or len(accession.strip()) < 4:
        return jsonify({
            'success': False,
            'error': 'Invalid accession ID'
        }), 400
    
    try:
        data = fetch_uniprot_protein_metadata(accession.strip())
        
        # Extract and normalize features
        features = []
        for feature in data.get('features', []):
            location = feature.get('location', {})
            start_data = location.get('start', {})
            end_data = location.get('end', {})
            
            # Extract numeric values from location objects
            start = start_data.get('value') if isinstance(start_data, dict) else start_data
            end = end_data.get('value') if isinstance(end_data, dict) else end_data
            
            if start and end:
                features.append({
                    'type': feature.get('type', 'Unknown'),
                    'description': feature.get('description', ''),
                    'location': {
                        'start': start,
                        'end': end
                    }
                })
        
        # Extract relevant information
        protein = {
            'accession': data.get('primaryAccession', accession),
            'name': (
                data.get('proteinDescription', {})
                .get('recommendedName', {})
                .get('fullName', {})
                .get('value')
                or data.get('proteinDescription', {})
                .get('submissionNames', [{}])[0]
                .get('fullName', {})
                .get('value')
                or 'Unknown protein'
            ),
            'organism': data.get('organism', {}).get('scientificName', 'Unknown organism'),
            'sequence': data.get('sequence', {}).get('value', ''),
            'length': data.get('sequence', {}).get('length', 0),
            'features': features
        }
        
        return jsonify({
            'success': True,
            'protein': protein
        }), 200
        
    except UniProtNotFound as e:
        return jsonify({
            'success': False,
            'error': f'Protein with accession "{accession}" not found in UniProt'
        }), 404
        
    except UniProtNetworkError as e:
        return jsonify({
            'success': False,
            'error': 'Failed to connect to UniProt'
        }), 503
        
    except Exception as e:
        return jsonify({
            'success': False,
            'error': 'Server error while fetching protein data'
        }), 500
