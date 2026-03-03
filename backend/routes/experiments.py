from flask import Blueprint, request, jsonify, session, send_file
import io
import matplotlib
matplotlib.use('Agg')  # headless backend — must be set before any pyplot import
from services.experiment_service import experiment_service
from services.experimental_data_parser import parser
from services.activity_calculator import activity_calculator
from services.sequence_analyzer import sequence_analyzer
from services.fingerprint_plot import resolve_structure, build_3d_fingerprint
from services.uniprot_client import fetch_uniprot_features_json
from models.experiment import Experiment, VariantData, Mutation, safe_float
from database import db
from sqlalchemy.orm import joinedload
import pandas as pd
import uuid
import math
import threading


def clean_dict_for_json(obj):
    """Convert NaN/Inf values to None for JSON serialization"""
    if isinstance(obj, dict):
        return {k: clean_dict_for_json(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [clean_dict_for_json(item) for item in obj]
    elif isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)):
        return None
    return obj


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

        # ── Validate FASTA structure before hitting UniProt ────────────────────
        # Give clear feedback for the two most common mistakes: pasting raw DNA
        # without a header, and uploading a non-DNA file.
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


@experiments_bp.route('/<experiment_id>/upload-data', methods=['POST'])
def upload_experimental_data(experiment_id: str):
    """
    Upload and process experimental data (TSV/JSON) for an experiment.
    Handles parsing, QC, sequence analysis, and activity score calculation.
    """
    user_id = require_auth()
    if not user_id:
        return jsonify({'success': False, 'error': 'Not authenticated'}), 401
    
    try:
        # Get experiment and verify ownership
        experiment = experiment_service.get_experiment_by_id(experiment_id, user_id)
        if not experiment:
            return jsonify({'success': False, 'error': 'Experiment not found'}), 404
        
        # Extract data we need from experiment (before closing session)
        exp_id = experiment.id
        wt_protein_seq = experiment.wt_protein_sequence
        plasmid_seq = experiment.plasmid_sequence
        
        # Get request data
        data = request.get_json()
        file_content = data.get('data', '')
        file_format = data.get('format', 'tsv').lower()
        
        if not file_content:
            return jsonify({'success': False, 'error': 'No file content provided'}), 400

        # ── Early content validation ──────────────────────────────────────────
        if file_format not in ('tsv', 'txt', 'json'):
            return jsonify({
                'success': False,
                'error': f'Unsupported file format "{file_format}". Please upload a .tsv or .json file.'
            }), 400

        if file_format in ('tsv', 'txt'):
            first_line = file_content.split('\n')[0] if file_content else ''
            if '\t' not in first_line:
                # Could be a CSV, Excel copy-paste, or wrong file entirely
                if ',' in first_line:
                    return jsonify({
                        'success': False,
                        'error': (
                            'The file appears to be comma-separated (CSV) rather than tab-separated (TSV). '
                            'Please export your data as a .tsv file with tab delimiters.'
                        )
                    }), 400
                else:
                    return jsonify({
                        'success': False,
                        'error': (
                            'The file does not appear to be tab-separated. '
                            'Please upload a .tsv file with columns separated by tab characters.'
                        )
                    }), 400

        if file_format == 'json':
            import json as _json
            try:
                _json.loads(file_content)
            except _json.JSONDecodeError as je:
                return jsonify({
                    'success': False,
                    'error': f'Invalid JSON file: {je.msg} at line {je.lineno}, column {je.colno}.'
                }), 400

        # Step 1: Parse and validate the data
        try:
            valid_df, control_df, rejected_df, parse_summary = parser.process_file(file_content, file_format)
            print(f"Parsed: {parse_summary['valid_rows']} valid, {parse_summary['control_rows']} controls, {parse_summary['rejected_rows']} rejected")
        except ValueError as e:
            print(f"Parse error: {e}")
            return jsonify({'success': False, 'error': str(e)}), 400

        if valid_df.empty and control_df.empty:
            return jsonify({
                'success': False,
                'error': 'No data rows could be parsed from the file.'
            }), 400
        
        # Step 2: Calculate activity scores (using both valid variants and controls)
        print("Step 2: Calculating activity scores...")
        try:
            # Combine valid variants with controls for activity calculation
            combined_df = pd.concat([valid_df, control_df], ignore_index=True) if not control_df.empty else valid_df
            # Calculate scores for the combined dataset
            scored_df = activity_calculator.calculate_activity_scores(combined_df)
            
            # Separate controls and variants after scoring
            control_scored_df = scored_df[scored_df['is_control'] == True].copy() if not control_df.empty else pd.DataFrame()
            valid_df = scored_df[scored_df['is_control'] == False].copy()
            
            print(f"Activity scores calculated for {len(valid_df)} variants and {len(control_scored_df)} controls")
        except ValueError as e:
            print(f"Activity calculation error: {e}")
            return jsonify({'success': False, 'error': f'Activity calculation failed: {str(e)}'}), 400
        
        # Step 3: Prepare records — sequence analysis runs separately via /analyze-sequences
        print("Step 3: Preparing records (sequence analysis deferred to /analyze-sequences)...")
        analyzed_variants = valid_df.to_dict('records')
        analyzed_controls = control_scored_df.to_dict('records') if not control_scored_df.empty else []

        
        # Step 4: Store in database in batches to avoid timeout
        # Store both variants AND controls so Gen 0 controls are available for mutation analysis
        print(f"Step 4: Storing {len(analyzed_variants)} variants and {len(analyzed_controls)} controls in database...")
        stored_count = 0
        metadata_columns = parse_summary['metadata_columns']
        BATCH_SIZE = 50  # Commit every 50 variants to avoid long transactions
        
        # Combine variants and controls for storage
        all_records = analyzed_variants + analyzed_controls
        
        for idx, variant_data in enumerate(all_records):
            # Create variant record
            variant = VariantData(
                experiment_id=exp_id,
                plasmid_variant_index=variant_data['plasmid_variant_index'],
                parent_plasmid_variant=variant_data.get('parent_plasmid_variant'),
                generation=int(variant_data['generation']),
                assembled_dna_sequence=variant_data['assembled_dna_sequence'],
                dna_yield=variant_data['dna_yield'],
                protein_yield=variant_data['protein_yield'],
                is_control=variant_data['is_control'],
                protein_sequence=variant_data.get('protein_sequence'),
                activity_score=variant_data.get('activity_score'),
                qc_status='passed',
                qc_message=None
            )
            
            # Store extra metadata
            if metadata_columns:
                metadata = {col: variant_data.get(col) for col in metadata_columns if col in variant_data}
                variant.extra_metadata = metadata
            
            # Store mutations using relationship (SQLAlchemy will handle variant_id)
            for mut in variant_data.get('mutations', []):
                mutation = Mutation(
                    position=mut['position'],
                    wild_type=mut['wt_aa'],
                    mutant=mut['mut_aa'],
                    mutation_type=mut.get('mutation_type', 'non-synonymous'),
                    generation_introduced=variant_data['generation']  # Simplified
                )
                variant.mutations.append(mutation)
            
            db.add(variant)
            stored_count += 1
            
            # Commit in batches to avoid long-running transactions
            if (idx + 1) % BATCH_SIZE == 0:
                db.commit()
                print(f"  Committed batch: {stored_count}/{len(all_records)} records")
        
        # Final commit for remaining items
        db.commit()
        print(f"Database commit successful. Stored {len(analyzed_variants)} variants + {len(analyzed_controls)} controls = {stored_count} total records.")
        
        # Step 5: Calculate statistics for response
        print("Step 5: Generating statistics...")
        generation_stats = activity_calculator.get_generation_statistics(valid_df)
        
        # Build response
        stats_list = generation_stats.to_dict('records') if not generation_stats.empty else []
        stats_list = clean_dict_for_json(stats_list)
        
        response = {
            'success': True,
            'parsed': parse_summary['total_rows'],
            'processed': stored_count,
            'variants': len(analyzed_variants),
            'controls': len(analyzed_controls),
            'passedQC': parse_summary['valid_rows'],
            'failedQC': parse_summary['rejected_rows'],
            'errors': [
                {
                    'row': r['qc_row_number'],
                    'message': r['qc_error_reason']
                }
                for r in parse_summary['rejected_details']
            ][:20],  # Limit to first 20 errors
            'warnings': [],
            'generationStats': stats_list
        }
        
        print(f"Upload completed successfully! Returning response.")
        return jsonify(response), 200
        
    except Exception as e:
        db.rollback()
        print(f"Error in upload_experimental_data: {str(e)}")
        import traceback
        traceback.print_exc()
        return jsonify({'success': False, 'error': 'Server error during data processing'}), 500


@experiments_bp.route('/<experiment_id>/variants', methods=['GET'])
def get_experiment_variants(experiment_id: str):
    """Get all variant data for an experiment"""
    user_id = require_auth()
    if not user_id:
        return jsonify({'success': False, 'error': 'Not authenticated'}), 401
    
    try:
        # Verify experiment ownership
        experiment = experiment_service.get_experiment_by_id(experiment_id, user_id)
        if not experiment:
            return jsonify({'success': False, 'error': 'Experiment not found'}), 404
        
        # Get variants with optional pagination
        exp_uuid = uuid.UUID(experiment_id) if isinstance(experiment_id, str) else experiment_id
        limit = min(int(request.args.get('limit', 1000)), 5000)  # Default 1000, max 5000
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
        
    except Exception as e:
        return jsonify({'success': False, 'error': 'Server error'}), 500


@experiments_bp.route('/<experiment_id>/top-performers', methods=['GET'])
def get_top_performers(experiment_id: str):
    """Get top performing variants by activity score"""
    user_id = require_auth()
    if not user_id:
        return jsonify({'success': False, 'error': 'Not authenticated'}), 401
    
    try:
        # Verify experiment ownership
        experiment = experiment_service.get_experiment_by_id(experiment_id, user_id)
        if not experiment:
            return jsonify({'success': False, 'error': 'Experiment not found'}), 404
        
        # Get top performers
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

        # ── Debug: log mutation counts so we can verify DB state ──────────────
        for v in variants:
            mut_list = v.mutations if include_mutations else []
            print(f"[TOP] variant pvi={v.plasmid_variant_index} gen={v.generation} "
                  f"mutations={len(mut_list)} activity={v.activity_score:.3f}")

        return jsonify({
            'success': True,
            'topPerformers': [v.to_dict(include_sequences=True, include_mutations=include_mutations) for v in variants]
        }), 200
        
    except Exception as e:
        return jsonify({'success': False, 'error': 'Server error'}), 500


@experiments_bp.route('/<experiment_id>/fingerprint/<variant_id>', methods=['GET'])
def get_mutation_fingerprint(experiment_id: str, variant_id: str):
    """
    Return the mutation fingerprint for a specific variant.
    Reconstructs the evolutionary lineage following parent_plasmid_variant
    backwards, then assigns generation_introduced for each mutation.
    This mirrors the teammate's mutation_fingerprinting.py algorithm but
    uses SQLAlchemy and the correct schema.
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

        # ── 1. Lightweight load of ALL variants (navigation only, no mutations) ─
        # plasmid_variant_index / parent_plasmid_variant / generation are all we
        # need to reconstruct the ancestry chain — this stays fast even with 300+
        # variants because we never touch the mutations relationship here.
        all_variants_light = (
            db.query(
                VariantData.id,
                VariantData.plasmid_variant_index,
                VariantData.parent_plasmid_variant,
                VariantData.generation,
                VariantData.activity_score,
            )
            .filter(VariantData.experiment_id == exp_uuid)
            .all()
        )

        # Build lookup tables (float keys are fine for PVI values like 0.0, 1.0…)
        by_pvi: dict = {}  # pvi -> row
        by_id:  dict = {}  # uuid -> row
        for row in all_variants_light:
            by_pvi[row.plasmid_variant_index] = row
            by_id[row.id] = row

        selected_light = by_id.get(var_uuid)
        if not selected_light:
            return jsonify({'success': False, 'error': 'Variant not found'}), 404

        # ── 2. Walk backwards from selected variant to root ────────────────────
        lineage_pvs = []   # list of (plasmid_variant_index, generation) in walk order
        current = selected_light
        visited: set = set()
        while current is not None and current.plasmid_variant_index not in visited:
            visited.add(current.plasmid_variant_index)
            lineage_pvs.append(current)
            parent_pvi = current.parent_plasmid_variant
            current = by_pvi.get(parent_pvi) if parent_pvi is not None and parent_pvi >= 0 else None

        # Sort chronologically (earliest generation first)
        lineage_pvs.sort(key=lambda r: r.generation)

        # ── 3. Load mutations ONLY for lineage variants ────────────────────────
        lineage_ids = [r.id for r in lineage_pvs]
        lineage_variants = (
            db.query(VariantData)
            .filter(VariantData.id.in_(lineage_ids))
            .options(joinedload(VariantData.mutations))
            .all()
        )
        lineage_by_id = {v.id: v for v in lineage_variants}
        # Preserve chronological order
        lineage = [lineage_by_id[r.id] for r in lineage_pvs if r.id in lineage_by_id]

        # ── 4. Delta walk: find which mutations are NEW at each generation ─────
        # "New" = present at this step but absent from the previous step.
        # This derives the correct generationIntroduced purely from the mutation
        # content at each generation — it does NOT rely on the stored
        # generation_introduced column (which was set incorrectly on the initial
        # upload and may not reflect true introduction generation).
        fingerprint: list = []
        seen_keys: set = set()       # avoid duplicating a mutation across branches
        prev_mut_keys: set = set()   # non-synonymous mut keys at previous step

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

        protein_length = len(experiment.wt_protein_sequence.strip()) if experiment.wt_protein_sequence else 0

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

    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({'success': False, 'error': 'Server error'}), 500


@experiments_bp.route('/<experiment_id>/fingerprint3d/<variant_id>', methods=['GET'])
def get_mutation_fingerprint_3d(experiment_id: str, variant_id: str):
    """
    Return a Plotly JSON figure for the circular mutation fingerprint with an
    optional 3D structure panel (AlphaFold / PDB via protein_accession).

    The figure JSON is consumed directly by react-plotly.js on the frontend.
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

        # ── 1. Load all variants (lightweight) for lineage walk ──────────────
        all_light = (
            db.query(
                VariantData.id,
                VariantData.plasmid_variant_index,
                VariantData.parent_plasmid_variant,
                VariantData.generation,
                VariantData.activity_score,
                VariantData.protein_sequence,
            )
            .filter(VariantData.experiment_id == exp_uuid)
            .all()
        )

        by_pvi = {r.plasmid_variant_index: r for r in all_light}
        by_id = {r.id: r for r in all_light}

        selected_light = by_id.get(var_uuid)
        if not selected_light:
            return jsonify({'success': False, 'error': 'Variant not found'}), 404

        # ── 2. Lineage walk ──────────────────────────────────────────────────
        lineage_pvs = []
        current = selected_light
        visited: set = set()
        while current is not None and current.plasmid_variant_index not in visited:
            visited.add(current.plasmid_variant_index)
            lineage_pvs.append(current)
            parent_pvi = current.parent_plasmid_variant
            current = by_pvi.get(parent_pvi) if parent_pvi is not None and parent_pvi >= 0 else None
        lineage_pvs.sort(key=lambda r: r.generation)

        # ── 3. Load mutations for the lineage ────────────────────────────────
        lineage_ids = [r.id for r in lineage_pvs]
        lineage_variants = (
            db.query(VariantData)
            .filter(VariantData.id.in_(lineage_ids))
            .options(joinedload(VariantData.mutations))
            .all()
        )
        lineage_by_id = {v.id: v for v in lineage_variants}
        lineage = [lineage_by_id[r.id] for r in lineage_pvs if r.id in lineage_by_id]

        # ── 4. Build fingerprint (delta walk) ────────────────────────────────
        # Collect ALL mutations (synonymous + non-synonymous) for the 3D view
        all_mutations_raw = {}
        prev_keys: set = set()

        for step in lineage:
            step_keys = {
                (m.position, m.wild_type, m.mutant, m.mutation_type)
                for m in (step.mutations or [])
            }
            new_keys = step_keys - {(k[0], k[1], k[2], k[3]) for k in prev_keys}
            for m in (step.mutations or []):
                key = (m.position, m.wild_type, m.mutant, m.mutation_type)
                if key in new_keys and key not in all_mutations_raw:
                    all_mutations_raw[key] = {
                        'position': m.position,
                        'wt_aa': m.wild_type,
                        'mut_aa': m.mutant,
                        'mutation_type': m.mutation_type,
                        'wt_codon': m.wt_codon or '',
                        'mut_codon': m.mut_codon or '',
                        'aa_change': f"{m.wild_type}{m.position}{m.mutant}",
                        'generation': step.generation,
                    }
            prev_keys = step_keys

        mutations = sorted(all_mutations_raw.values(), key=lambda m: m['position'])

        # ── 5. Mutation counts ─────────────────────────────────────────────
        num_synonymous = sum(1 for m in mutations if m['mutation_type'] == 'synonymous')
        num_nonsynonymous = len(mutations) - num_synonymous

        # ── 6. Resolve 3D structure ──────────────────────────────────────────
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

        # ── 7. UniProt feature annotations ───────────────────────────────────
        feature_annotations: list = []
        if uniprot_id:
            try:
                feature_annotations = fetch_uniprot_features_json(uniprot_id)
            except Exception:
                pass

        # ── 8. Build Plotly figure and return as JSON ────────────────────────
        variant_pvi = selected_light.plasmid_variant_index
        fig = build_3d_fingerprint(
            lineage_mutations=mutations,
            backbone=backbone,
            pdb_coords=pdb_coords,
            variant_id=variant_pvi,
            uniprot_id=uniprot_id,
            structure_source=structure_info.get('source', structure_info.get('status', 'Structure')),
            feature_annotations=feature_annotations,
        )
        import json as _json
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

    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({'success': False, 'error': 'Server error'}), 500


def get_experiment_statistics(experiment_id: str):
    """Get statistical summary of experiment data"""
    user_id = require_auth()
    if not user_id:
        return jsonify({'success': False, 'error': 'Not authenticated'}), 401
    
    try:
        # Verify experiment ownership
        experiment = experiment_service.get_experiment_by_id(experiment_id, user_id)
        if not experiment:
            return jsonify({'success': False, 'error': 'Experiment not found'}), 404
        
        # Get all variants for this experiment
        exp_uuid = uuid.UUID(experiment_id) if isinstance(experiment_id, str) else experiment_id
        variants = db.query(VariantData).filter_by(
            experiment_id=exp_uuid
        ).all()
        
        if not variants:
            return jsonify({
                'success': True,
                'statistics': {
                    'totalVariants': 0,
                    'passedQC': 0,
                    'failedQC': 0,
                    'generations': []
                }
            }), 200
        
        # Convert to DataFrame for analysis
        variants_data = []
        for v in variants:
            variants_data.append({
                'generation': v.generation,
                'activity_score': v.activity_score,
                'is_control': v.is_control,
                'qc_status': v.qc_status
            })
        
        df = pd.DataFrame(variants_data)
        
        # Calculate generation statistics
        generation_stats = activity_calculator.get_generation_statistics(df)
        stats_list = generation_stats.to_dict('records') if not generation_stats.empty else []
        stats_list = clean_dict_for_json(stats_list)
        
        return jsonify({
            'success': True,
            'statistics': {
                'totalVariants': len(variants),
                'passedQC': len([v for v in variants if v.qc_status == 'passed']),
                'failedQC': len([v for v in variants if v.qc_status == 'failed']),
                'generationStats': stats_list
            }
        }), 200
        
    except Exception as e:
        print(f"Error in get_experiment_statistics: {str(e)}")
        return jsonify({'success': False, 'error': 'Server error'}), 500


def _set_analysis_status(exp_id, status: str, message: str):
    """Helper to update experiment analysis status safely."""
    try:
        experiment = db.query(Experiment).filter_by(id=exp_id).first()
        if experiment:
            experiment.analysis_status = status
            experiment.analysis_message = message
            db.commit()
    except Exception as e:
        print(f"Could not update analysis status: {e}")
        db.rollback()


def _run_analysis_background(app, exp_id, wt_protein_seq: str, plasmid_seq: str):
    """
    Runs in a background thread. Uses app context so the DB session works correctly.
    The HTTP response has already been returned before this runs.
    """
    with app.app_context():
        try:
            print(f"[BG] Starting sequence analysis for {exp_id}...")

            all_data = db.query(VariantData).filter_by(experiment_id=exp_id).all()
            if not all_data:
                _set_analysis_status(exp_id, 'failed', 'No data found to analyze')
                return

            controls = [v for v in all_data if v.is_control and v.generation == 0]
            variants = [v for v in all_data if not v.is_control]
            print(f"[BG] {len(variants)} variants, {len(controls)} Gen-0 controls")

            # Use Gen-0 control plasmid as WT reference if available
            if controls:
                print("[BG] Using Generation 0 control as WT reference...")
                ref_plasmid = controls[0].assembled_dna_sequence
            else:
                print("[BG] No Gen-0 controls found — using experiment plasmid as reference")
                ref_plasmid = plasmid_seq

            variants_data = [
                {'id': str(v.id), 'assembled_dna_sequence': v.assembled_dna_sequence, 'generation': v.generation}
                for v in variants
            ]

            print("[BG] Running sequence analyzer...")
            analyzed = sequence_analyzer.analyze_variant_batch(
                variants_data, wt_protein_seq, ref_plasmid
            )
            print(f"[BG] Analysis complete for {len(analyzed)} variants")

            # Build dicts from already-loaded ORM objects — avoids 1 query per variant
            variant_dict = {str(v.id): v for v in variants}
            variant_ids = [v.id for v in variants]

            # ── Correct generation_introduced via lineage walk ─────────────────
            # For each variant, mutations that appear in its parent are *inherited*
            # and keep the parent's generation. Mutations not in the parent are new.
            #
            # Build: { plasmid_variant_index (float) -> { mut_key -> generation } }
            # Process in ascending generation order.
            idx_to_mutations: dict = {}   # pvi -> {(position, wt, mut): generation}
            pvi_to_variant   = {v.plasmid_variant_index: v for v in variants}

            analyzed_sorted = sorted(analyzed, key=lambda r: variant_dict[r['id']].generation)

            for result in analyzed_sorted:
                v = variant_dict[result['id']]
                parent_pvi  = v.parent_plasmid_variant
                parent_muts = idx_to_mutations.get(parent_pvi, {})

                this_muts: dict = {}
                for mut in result.get('mutations', []):
                    key = (mut['position'], mut['wt_aa'], mut['mut_aa'])
                    if key in parent_muts:
                        # Inherited — keep the generation it was first introduced
                        this_muts[key] = parent_muts[key]
                    else:
                        # New this generation
                        this_muts[key] = v.generation
                idx_to_mutations[v.plasmid_variant_index] = this_muts

            # ── Build new Mutation objects and protein updates (no DB writes yet) ────
            # ALL computation completes before ANY destructive DB operation.
            # This means a server restart or exception during analysis leaves the
            # existing mutations intact — the old delete-first approach wiped them.
            new_mutations = []
            protein_updates = {}   # variant_id -> new protein_sequence string
            BATCH_SIZE = 200

            for result in analyzed:
                variant = variant_dict.get(result['id'])
                if not variant:
                    continue

                protein_updates[result['id']] = result.get('protein_sequence')

                gen_map = idx_to_mutations.get(variant.plasmid_variant_index, {})
                for mut in result.get('mutations', []):
                    key = (mut['position'], mut['wt_aa'], mut['mut_aa'])
                    gen_introduced = gen_map.get(key, variant.generation)
                    new_mutations.append(Mutation(
                        variant_id=variant.id,
                        position=mut['position'],
                        wild_type=mut['wt_aa'],
                        mutant=mut['mut_aa'],
                        wt_codon=mut.get('wt_codon'),
                        mut_codon=mut.get('mut_codon'),
                        mut_aa=mut.get('mut_aa'),
                        mutation_type=mut.get('mutation_type', 'non-synonymous'),
                        generation_introduced=gen_introduced
                    ))

            print(f"[BG] Built {len(new_mutations)} mutations for {len(protein_updates)} variants — writing to DB...")

            # ── Atomic swap: delete old mutations then insert new ones ─────────
            # Both steps are in the same transaction so a crash leaves the DB in
            # a consistent state (either fully old or fully new data).
            updated_count = 0
            if variant_ids:
                deleted = db.query(Mutation).filter(
                    Mutation.variant_id.in_(variant_ids)
                ).delete(synchronize_session=False)
                print(f"[BG] Deleted {deleted} old mutations")

            for result_id, new_protein in protein_updates.items():
                variant = variant_dict.get(result_id)
                if variant:
                    variant.protein_sequence = new_protein
                    updated_count += 1
                    # Flush protein-sequence updates in batches to keep memory low
                    if updated_count % BATCH_SIZE == 0:
                        db.flush()
                        print(f"[BG] Flushed protein sequences: {updated_count}/{len(protein_updates)}")

            if new_mutations:
                # add_all() is used instead of bulk_save_objects() — it respects
                # the current transaction and fires ORM events correctly.
                db.add_all(new_mutations)
                print(f"[BG] Staged {len(new_mutations)} mutations for insert")

            db.commit()
            print(f"[BG] Database update complete: {updated_count} variants, {len(new_mutations)} mutations")
            _set_analysis_status(exp_id, 'completed', f'Successfully analyzed {updated_count} variants')

        except Exception as e:
            import traceback
            traceback.print_exc()
            db.rollback()
            _set_analysis_status(exp_id, 'failed', f'Analysis failed: {str(e)}')
        finally:
            # Always release the scoped session for this thread back to the pool.
            # Without this, the background-thread session stays in the registry
            # and holds a DB connection open indefinitely.
            db.remove()


@experiments_bp.route('/<experiment_id>/analyze-sequences', methods=['POST'])
def analyze_sequences(experiment_id: str):
    """
    Kick off background sequence analysis. Returns 200 immediately.
    Poll GET /<experiment_id> to check analysis_status field.
    """
    user_id = require_auth()
    if not user_id:
        return jsonify({'success': False, 'error': 'Not authenticated'}), 401

    try:
        experiment = experiment_service.get_experiment_by_id(experiment_id, user_id)
        if not experiment:
            return jsonify({'success': False, 'error': 'Experiment not found'}), 404

        if experiment.analysis_status == 'analyzing':
            return jsonify({'success': False, 'error': 'Analysis already in progress'}), 409

        exp_id = experiment.id
        wt_protein_seq = experiment.wt_protein_sequence
        plasmid_seq = experiment.plasmid_sequence

        # Set status synchronously before spawning thread
        experiment.analysis_status = 'analyzing'
        experiment.analysis_message = 'Analysis queued...'
        db.commit()

        # Spawn background thread — HTTP response returns immediately
        from flask import current_app
        app = current_app._get_current_object()

        thread = threading.Thread(
            target=_run_analysis_background,
            args=(app, exp_id, wt_protein_seq, plasmid_seq),
            daemon=True
        )
        thread.start()

        return jsonify({
            'success': True,
            'message': 'Sequence analysis started',
            'status': 'analyzing'
        }), 200

    except Exception as e:
        import traceback
        traceback.print_exc()
        db.rollback()
        return jsonify({'success': False, 'error': f'Failed to start analysis: {str(e)}'}), 500


@experiments_bp.route('/<experiment_id>/plots/activity-distribution', methods=['GET'])
def plot_activity_distribution(experiment_id: str):
    """Return a PNG of the activity score distribution violin plot.
    Runs the original matplotlib/seaborn code from activity_score_per_gen.py
    server-side and streams the image back.
    """
    user_id = session.get('user_id')
    if not user_id:
        return jsonify({'success': False, 'error': 'Unauthorized'}), 401

    try:
        import matplotlib.pyplot as plt
        import numpy as np

        experiment = experiment_service.get_experiment_by_id(experiment_id, user_id)
        if not experiment:
            return jsonify({'success': False, 'error': 'Experiment not found'}), 404

        exp_uuid = uuid.UUID(experiment_id) if isinstance(experiment_id, str) else experiment_id
        variants = (
            db.query(VariantData)
            .filter(
                VariantData.experiment_id == exp_uuid,
                VariantData.qc_status == 'passed',
                VariantData.activity_score.isnot(None),
            )
            .all()
        )

        if not variants:
            return jsonify({'success': False, 'error': 'No QC-passed variants with activity scores'}), 404

        df = pd.DataFrame([
            {
                'Directed_Evolution_Generation': v.generation,
                'activity_score': float(v.activity_score),
            }
            for v in variants
        ])
        df = df.dropna(subset=['Directed_Evolution_Generation', 'activity_score'])

        # ── original activity_score_per_gen.py plot logic ─────────────────────
        gens = sorted(df['Directed_Evolution_Generation'].unique())

        # Filter out generations with fewer than 2 points (violinplot requires ≥2)
        data = [
            df[df['Directed_Evolution_Generation'] == g]['activity_score'].values.astype(float)
            for g in gens
            if len(df[df['Directed_Evolution_Generation'] == g]) >= 2
        ]
        gens = [
            g for g in gens
            if len(df[df['Directed_Evolution_Generation'] == g]) >= 2
        ]

        if not data:
            return jsonify({'success': False, 'error': 'Not enough data points per generation to draw violin plots (need ≥2 per generation)'}), 404

        fig, ax = plt.subplots(figsize=(10, 7))
        fig.patch.set_facecolor('#0f172a')   # dark background to match portal
        ax.set_facecolor('#1e293b')

        violins = ax.violinplot(data, widths=0.9, points=200)

        for pc in violins['bodies']:
            pc.set_facecolor('#add8e6')
            pc.set_edgecolor('#5aa0c8')
            pc.set_alpha(0.8)
        for part in ('cbars', 'cmins', 'cmaxes'):
            if part in violins:
                violins[part].set_visible(False)

        for i, v in enumerate(data, start=1):
            vmin = np.min(v)
            vmax = np.max(v)
            mean = np.mean(v)
            median = np.median(v)

            ax.vlines(i, vmin, vmax, linewidth=2, colors='#1e40af')
            ax.hlines(vmin,   i - 0.12, i + 0.12, linewidth=2, colors='#1e40af')
            ax.hlines(vmax,   i - 0.12, i + 0.12, linewidth=2, colors='#1e40af')
            ax.hlines(mean,   i - 0.18, i + 0.18, linewidth=3, colors='#1e40af')
            ax.hlines(median, i - 0.18, i + 0.18, linewidth=3, colors='#1e40af')

        ax.set_xticks(range(1, len(gens) + 1))
        ax.set_xticklabels([f'Gen {int(g)}' for g in gens], color='#94a3b8')
        ax.tick_params(colors='#94a3b8')
        ax.set_xlabel('Generation', color='#94a3b8')
        ax.set_ylabel('Activity Score', color='#94a3b8')
        ax.set_title('Activity Score Distribution by Generation', color='#e2e8f0')
        for spine in ax.spines.values():
            spine.set_color('#334155')
        ax.grid(True, color='#334155', linewidth=0.5)

        plt.tight_layout()
        # ───────────────────────────────────────────────────────────────────────

        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=150, bbox_inches='tight',
                    facecolor=fig.get_facecolor())
        plt.close(fig)
        buf.seek(0)

        return send_file(buf, mimetype='image/png')

    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({'success': False, 'error': str(e)}), 500
