"""
Data upload routes.

Routes registered here:
  POST  /api/experiments/<experiment_id>/preview-mapping
  POST  /api/experiments/<experiment_id>/upload-data
"""
from flask import request, jsonify
import pandas as pd
import math as _math

from database import db
from models.experiment import VariantData, Mutation
from services.experiment_service import experiment_service
from services.experimental_data_parser import parser
from services.activity_calculator import activity_calculator
from ._base import experiments_bp, require_auth, clean_dict_for_json


@experiments_bp.route('/<experiment_id>/preview-mapping', methods=['POST'])
def preview_column_mapping(experiment_id: str):
    """
    Preview the auto-detected column mapping for an uploaded file without
    persisting any data.  Used by the frontend mapping-confirmation step.

    Request body (JSON):
      data    (str) — raw file text
      format  (str) — "tsv" or "json"

    Response 200:
      { success, raw_columns, column_mapping, metadata_columns,
        missing_required, canonical_fields }
    """
    user_id = require_auth()
    if not user_id:
        return jsonify({'success': False, 'error': 'Not authenticated'}), 401

    body = request.get_json(silent=True) or {}
    file_content: str = body.get('data', '')
    file_format: str = body.get('format', 'tsv').lower()

    if not file_content:
        return jsonify({'success': False, 'error': 'No file content provided'}), 400

    try:
        preview = parser.preview_mapping(file_content, file_format)
        return jsonify({'success': True, **preview}), 200
    except ValueError as exc:
        return jsonify({'success': False, 'error': str(exc)}), 400
    except Exception:
        return jsonify({'success': False, 'error': 'Failed to parse file header'}), 500


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
        # Optional user-confirmed mapping from the frontend review step
        column_mapping_override = data.get('column_mapping') or None

        if not file_content:
            return jsonify({'success': False, 'error': 'No file content provided'}), 400

        # ── Early content validation ────────────────────────────────────────
        if file_format not in ('tsv', 'txt', 'json'):
            return jsonify({
                'success': False,
                'error': f'Unsupported file format "{file_format}". Please upload a .tsv or .json file.'
            }), 400

        if file_format in ('tsv', 'txt'):
            first_line = file_content.split('\n')[0] if file_content else ''
            if '\t' not in first_line:
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
            valid_df, control_df, rejected_df, parse_summary = parser.process_file(
                file_content, file_format,
                column_mapping_override=column_mapping_override,
            )
            print(f"Parsed: {parse_summary['valid_rows']} valid, "
                  f"{parse_summary['control_rows']} controls, "
                  f"{parse_summary['rejected_rows']} rejected")

            combined_check = pd.concat([valid_df, control_df], ignore_index=True) if not control_df.empty else valid_df
            if combined_check.empty:
                raise ValueError(
                    "Uploaded file contains no data rows. "
                    "Provide at least one experimental record."
                )

        except ValueError as e:
            print(f"Parse error: {e}")
            return jsonify({'success': False, 'error': str(e)}), 400

        if control_df.empty:
            return jsonify({
                'success': False,
                'error': (
                    'No control rows were found in the uploaded data. '
                    'At least one row must be marked as a control '
                    '(is_control = TRUE / 1) before activity scores can be calculated.'
                )
            }), 400

        # Step 2: Calculate activity scores
        print("Step 2: Calculating activity scores...")
        try:
            combined_df = (
                pd.concat([valid_df, control_df], ignore_index=True)
                if not control_df.empty else valid_df
            )
            scored_df = activity_calculator.calculate_activity_scores(combined_df)

            control_scored_df = (
                scored_df[scored_df['is_control'] == True].copy()
                if not control_df.empty else pd.DataFrame()
            )
            valid_df = scored_df[scored_df['is_control'] == False].copy()

            print(f"Activity scores calculated for {len(valid_df)} variants "
                  f"and {len(control_scored_df)} controls")
        except ValueError as e:
            print(f"Activity calculation error: {e}")
            return jsonify({'success': False, 'error': f'Activity calculation failed: {str(e)}'}), 400

        # Step 3: Prepare records (sequence analysis deferred to /analyze-sequences)
        print("Step 3: Preparing records...")
        analyzed_variants = valid_df.to_dict('records')
        analyzed_controls = control_scored_df.to_dict('records') if not control_scored_df.empty else []

        # Step 4: Store in database in batches
        print(f"Step 4: Storing {len(analyzed_variants)} variants and "
              f"{len(analyzed_controls)} controls in database...")
        stored_count = 0
        metadata_columns = parse_summary['metadata_columns']
        BATCH_SIZE = 50

        all_records = analyzed_variants + analyzed_controls

        for idx, variant_data in enumerate(all_records):
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

            # Store extra metadata — coerce NaN/inf to None for JSON safety
            if metadata_columns:
                metadata = {}
                for col in metadata_columns:
                    if col not in variant_data:
                        continue
                    val = variant_data[col]
                    if isinstance(val, float) and (_math.isnan(val) or _math.isinf(val)):
                        val = None
                    metadata[col] = val
                variant.extra_metadata = metadata or None

            for mut in variant_data.get('mutations', []):
                mutation = Mutation(
                    position=mut['position'],
                    wild_type=mut['wt_aa'],
                    mutant=mut['mut_aa'],
                    mutation_type=mut.get('mutation_type', 'non-synonymous'),
                    generation_introduced=variant_data['generation']
                )
                variant.mutations.append(mutation)

            db.add(variant)
            stored_count += 1

            if (idx + 1) % BATCH_SIZE == 0:
                db.commit()
                print(f"  Committed batch: {stored_count}/{len(all_records)} records")

        db.commit()
        print(f"Database commit successful. "
              f"Stored {len(analyzed_variants)} variants + {len(analyzed_controls)} controls "
              f"= {stored_count} total records.")

        # Step 5: Calculate statistics for response
        print("Step 5: Generating statistics...")
        generation_stats = activity_calculator.get_generation_statistics(valid_df)

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
                {'row': r['qc_row_number'], 'message': r['qc_error_reason']}
                for r in parse_summary['rejected_details']
            ][:20],
            'warnings': [],
            'generationStats': stats_list
        }

        print("Upload completed successfully!")
        return jsonify(response), 200

    except Exception as e:
        db.rollback()
        print(f"Error in upload_experimental_data: {str(e)}")
        import traceback
        traceback.print_exc()
        return jsonify({'success': False, 'error': 'Server error during data processing'}), 500
