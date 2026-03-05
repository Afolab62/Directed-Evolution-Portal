from flask import Blueprint, jsonify, Response
from services.uniprot_client import (
    fetch_uniprot_protein_metadata,
    fetch_uniprot_detailed,
    fetch_uniprot_fasta,
    UniProtNotFound,
    UniProtNetworkError,
    UniProtError,
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
        detail = fetch_uniprot_detailed(accession.strip())

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

        # Extract relevant information, merging in the enriched detail fields
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
            'organism':          data.get('organism', {}).get('scientificName', 'Unknown organism'),
            'sequence':          data.get('sequence', {}).get('value', ''),
            'length':            data.get('sequence', {}).get('length', 0),
            'features':          features,
            # Enriched fields from detailed fetch
            'gene_names':        detail['gene_names'],
            'keywords':          detail['keywords'],
            'taxonomy_lineage':  detail['taxonomy_lineage'],
            'pdb_ids':           detail['pdb_ids'],
            'alphafold_id':      detail['alphafold_id'],
            'interpro_ids':      detail['interpro_ids'],
            'pfam_ids':          detail['pfam_ids'],
            'kegg_ids':          detail['kegg_ids'],
            'go_terms':          detail['go_terms'],
            'crossrefs':         detail['crossrefs'],
        }

        return jsonify({
            'success': True,
            'protein': protein
        }), 200

    except UniProtNotFound:
        return jsonify({
            'success': False,
            'error': f'Protein with accession "{accession}" not found in UniProt'
        }), 404

    except UniProtNetworkError:
        return jsonify({
            'success': False,
            'error': 'Failed to connect to UniProt'
        }), 503

    except Exception:
        import traceback
        traceback.print_exc()
        return jsonify({
            'success': False,
            'error': 'Server error while fetching protein data'
        }), 500


@uniprot_bp.route('/<accession>/fasta', methods=['GET'])
def download_fasta(accession: str):
    """
    Download the raw UniProt FASTA for an accession as a plain-text file.

    The browser will be prompted to save it as ``<accession>.fasta``.
    """
    if not accession or len(accession.strip()) < 4:
        return jsonify({'success': False, 'error': 'Invalid accession ID'}), 400

    try:
        fasta_text = fetch_uniprot_fasta(accession.strip())
        return Response(
            fasta_text,
            status=200,
            mimetype='text/plain',
            headers={
                'Content-Disposition': f'attachment; filename="{accession.strip()}.fasta"',
                'Content-Type': 'text/plain; charset=utf-8',
            },
        )

    except UniProtNotFound:
        return jsonify({
            'success': False,
            'error': f'Protein with accession "{accession}" not found in UniProt'
        }), 404

    except UniProtNetworkError:
        return jsonify({
            'success': False,
            'error': 'Failed to connect to UniProt — check network connectivity'
        }), 503

    except UniProtError as e:
        return jsonify({'success': False, 'error': str(e)}), 502

    except Exception:
        return jsonify({'success': False, 'error': 'Server error while downloading FASTA'}), 500
