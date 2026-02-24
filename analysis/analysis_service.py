"""
app/services/analysis_service.py
---------------------------------
Orchestration layer — the ONLY place where database queries and
analysis classes meet.

Responsibility
--------------
1. Fetch required data from the DB (via database/queries.py)
2. Construct SequenceAnalyzer with the injected StagingResult
3. Run analysis for each variant
4. Return results to the calling Flask route

Neither SequenceAnalyzer nor database/queries.py know about each other.
This service is what connects them.
"""

import logging
import sqlite3
from typing import Optional

from database.queries import (
    get_staging_result_by_experiment,
    get_variant_plasmid_sequence,
    get_wt_plasmid_sequence,
    get_wt_protein_sequence,
)
from analysis.sequence_analyzer import SequenceAnalyzer
from database.models import StagingResult

logger = logging.getLogger(__name__)


def run_variant_analysis(
    experiment_id: int,
    plasmid_variant_index: str,
    db_path: str = "database.db",
    wt_protein_seq: Optional[str] = None,
    wt_plasmid_seq: Optional[str] = None,
    variant_plasmid_seq: Optional[str] = None,
) -> dict:
    """
    Full analysis pipeline for a single variant.

    Sequences can be passed in directly (e.g. during parsing when they
    are already in memory) OR left as None and fetched from the DB.
    Passing them in avoids redundant DB round-trips when you already
    have the data available.

    Parameters
    ----------
    experiment_id : int
        ID of the parent experiment (used to fetch staging result + WT seqs).
    plasmid_variant_index : str
        Unique identifier of the variant being analysed.
    db_path : str
        Path to the SQLite database file.
    wt_protein_seq : str, optional
        WT protein sequence. Fetched from DB if not provided.
    wt_plasmid_seq : str, optional
        WT plasmid sequence. Fetched from DB if not provided.
    variant_plasmid_seq : str, optional
        Variant plasmid sequence. Fetched from DB if not provided.

    Returns
    -------
    dict with analysis results from SequenceAnalyzer.analyze_variant()
    plus 'plasmid_variant_index' and 'staging_result_id' for traceability.

    Raises
    ------
    ValueError  if required sequences are missing or staging result is invalid.
    """

    # ------------------------------------------------------------------ #
    # Step 1 — Fetch staging result (always from DB, it's always needed)  #
    # ------------------------------------------------------------------ #
    staging_result: Optional[StagingResult] = get_staging_result_by_experiment(
        experiment_id=experiment_id,
        db_path=db_path,
    )

    if staging_result is None:
        logger.warning(
            f"No valid staging result for experiment_id={experiment_id}. "
            "SequenceAnalyzer will fall back to Smith-Waterman — this is slower."
        )

    # ------------------------------------------------------------------ #
    # Step 2 — Resolve sequences (use provided values or fetch from DB)   #
    # ------------------------------------------------------------------ #
    if wt_protein_seq is None:
        wt_protein_seq = get_wt_protein_sequence(experiment_id, db_path)
        if wt_protein_seq is None:
            raise ValueError(
                f"WT protein sequence not found for experiment_id={experiment_id}. "
                "Ensure staging completed successfully and the UniProt sequence is stored."
            )

    if wt_plasmid_seq is None:
        wt_plasmid_seq = get_wt_plasmid_sequence(experiment_id, db_path)
        if wt_plasmid_seq is None:
            raise ValueError(
                f"WT plasmid sequence not found for experiment_id={experiment_id}."
            )

    if variant_plasmid_seq is None:
        variant_plasmid_seq = get_variant_plasmid_sequence(plasmid_variant_index, db_path)
        if variant_plasmid_seq is None:
            raise ValueError(
                f"Variant plasmid sequence not found for "
                f"Plasmid_Variant_Index={plasmid_variant_index}."
            )

    # ------------------------------------------------------------------ #
    # Step 3 — Build analyser (staging result injected, no DB calls here) #
    # ------------------------------------------------------------------ #
    analyzer = SequenceAnalyzer(
        wt_protein_seq=wt_protein_seq,
        wt_plasmid_seq=wt_plasmid_seq,
        staging_result=staging_result,   # None triggers SW fallback automatically
    )

    # ------------------------------------------------------------------ #
    # Step 4 — Run analysis and return                                    #
    # ------------------------------------------------------------------ #
    result = analyzer.analyze_variant(variant_plasmid_seq)

    result['plasmid_variant_index'] = plasmid_variant_index
    result['staging_result_id'] = staging_result.id if staging_result else None

    logger.info(
        f"Analysis complete for variant={plasmid_variant_index}: "
        f"{result['num_mutations']} mutations "
        f"({result['num_synonymous']} syn, {result['num_nonsynonymous']} non-syn)"
    )

    return result


def run_batch_analysis(
    experiment_id: int,
    variant_indices: list,
    db_path: str = "database.db",
    wt_protein_seq: Optional[str] = None,
    wt_plasmid_seq: Optional[str] = None,
) -> list:
    """
    Run analysis across a batch of variants for an experiment.

    Fetches the staging result and WT sequences ONCE, then reuses the same
    SequenceAnalyzer instance for all variants — much faster than calling
    run_variant_analysis() in a loop.

    Parameters
    ----------
    experiment_id : int
    variant_indices : list of Plasmid_Variant_Index strings
    db_path : str
    wt_protein_seq : str, optional — fetched from DB if not provided
    wt_plasmid_seq : str, optional — fetched from DB if not provided

    Returns
    -------
    List of result dicts (same format as run_variant_analysis).
    Failed variants are included with an 'error' key instead of sequence data.
    """

    # Fetch shared resources once
    staging_result = get_staging_result_by_experiment(experiment_id, db_path)

    if wt_protein_seq is None:
        wt_protein_seq = get_wt_protein_sequence(experiment_id, db_path)
    if wt_plasmid_seq is None:
        wt_plasmid_seq = get_wt_plasmid_sequence(experiment_id, db_path)

    if not wt_protein_seq or not wt_plasmid_seq:
        raise ValueError(
            f"Missing WT sequences for experiment_id={experiment_id}. "
            "Cannot proceed with batch analysis."
        )

    # Build a single shared analyser instance
    analyzer = SequenceAnalyzer(
        wt_protein_seq=wt_protein_seq,
        wt_plasmid_seq=wt_plasmid_seq,
        staging_result=staging_result,
    )

    results = []

    for variant_index in variant_indices:
        try:
            variant_seq = get_variant_plasmid_sequence(variant_index, db_path)

            if variant_seq is None:
                raise ValueError(f"Sequence not found for variant={variant_index}.")

            result = analyzer.analyze_variant(variant_seq)
            result['plasmid_variant_index'] = variant_index
            result['staging_result_id'] = staging_result.id if staging_result else None
            results.append(result)

        except Exception as exc:
            logger.error(f"Analysis failed for variant={variant_index}: {exc}")
            results.append({
                'plasmid_variant_index': variant_index,
                'error': str(exc),
            })

    logger.info(
        f"Batch analysis complete: {len(results)} variants processed "
        f"for experiment_id={experiment_id}."
    )

    return results
