"""
Sequence analysis routes and background worker.

Routes registered here:
  POST  /api/experiments/<experiment_id>/analyze-sequences
"""
import threading

from flask import request, jsonify, current_app

from database import db
from models.experiment import Experiment, VariantData, Mutation
from services.experiment_service import experiment_service
from services.sequence_analyzer import sequence_analyzer
from ._base import experiments_bp, require_auth


# ──────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ──────────────────────────────────────────────────────────────────────────────

def _set_analysis_status(exp_id, status: str, message: str):
    """Update experiment analysis_status safely (used from background thread)."""
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
    Run sequence analysis in a background thread.
    Uses an explicit app context so the DB session works correctly.
    The HTTP response has already been returned before this executes.
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
                {
                    'id': str(v.id),
                    'assembled_dna_sequence': v.assembled_dna_sequence,
                    'generation': v.generation,
                }
                for v in variants
            ]

            print("[BG] Running sequence analyzer...")
            analyzed = sequence_analyzer.analyze_variant_batch(
                variants_data, wt_protein_seq, ref_plasmid
            )
            print(f"[BG] Analysis complete for {len(analyzed)} variants")

            variant_dict = {str(v.id): v for v in variants}
            variant_ids  = [v.id for v in variants]

            # ── Derive generation_introduced via lineage walk ─────────────────
            # { plasmid_variant_index -> { (pos, wt, mut): generation } }
            idx_to_mutations: dict = {}

            analyzed_sorted = sorted(
                analyzed, key=lambda r: variant_dict[r['id']].generation
            )

            for result in analyzed_sorted:
                v = variant_dict[result['id']]
                parent_pvi  = v.parent_plasmid_variant
                parent_muts = idx_to_mutations.get(parent_pvi, {})

                this_muts: dict = {}
                for mut in result.get('mutations', []):
                    key = (mut['position'], mut['wt_aa'], mut['mut_aa'])
                    this_muts[key] = parent_muts.get(key, v.generation)
                idx_to_mutations[v.plasmid_variant_index] = this_muts

            # ── Build new Mutation objects (all computation before any DB write) ─
            new_mutations = []
            protein_updates = {}
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
                        generation_introduced=gen_introduced,
                    ))

            print(f"[BG] Built {len(new_mutations)} mutations for "
                  f"{len(protein_updates)} variants — writing to DB...")

            # ── Atomic replace: delete old → insert new ────────────────────────
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
                    if updated_count % BATCH_SIZE == 0:
                        db.flush()
                        print(f"[BG] Flushed protein sequences: "
                              f"{updated_count}/{len(protein_updates)}")

            if new_mutations:
                db.add_all(new_mutations)
                print(f"[BG] Staged {len(new_mutations)} mutations for insert")

            db.commit()
            print(f"[BG] Database update complete: "
                  f"{updated_count} variants, {len(new_mutations)} mutations")
            _set_analysis_status(
                exp_id, 'completed',
                f'Successfully analyzed {updated_count} variants'
            )

        except Exception as e:
            import traceback
            traceback.print_exc()
            db.rollback()
            _set_analysis_status(exp_id, 'failed', f'Analysis failed: {str(e)}')
        finally:
            # Release the scoped session for this thread back to the pool.
            db.remove()


# ──────────────────────────────────────────────────────────────────────────────
# Route
# ──────────────────────────────────────────────────────────────────────────────

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

        exp_id         = experiment.id
        wt_protein_seq = experiment.wt_protein_sequence
        plasmid_seq    = experiment.plasmid_sequence

        # Update status synchronously before spawning thread
        experiment.analysis_status  = 'analyzing'
        experiment.analysis_message = 'Analysis queued...'
        db.commit()

        app = current_app._get_current_object()
        thread = threading.Thread(
            target=_run_analysis_background,
            args=(app, exp_id, wt_protein_seq, plasmid_seq),
            daemon=True,
        )
        thread.start()

        return jsonify({
            'success': True,
            'message': 'Sequence analysis started',
            'status': 'analyzing',
        }), 200

    except Exception as e:
        import traceback
        traceback.print_exc()
        db.rollback()
        return jsonify({'success': False, 'error': f'Failed to start analysis: {str(e)}'}), 500
