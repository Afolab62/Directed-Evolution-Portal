"""
database/queries.py
-------------------
All raw SQL queries for the project live here.
No analysis logic — just fetch, return typed objects.

This is the ONLY file that should import sqlite3 and talk to the DB
for read operations. Write operations (INSERT/UPDATE) belong in
database/write.py (to be added as the project grows).
"""

import sqlite3
import logging
from typing import Optional

from database.models import StagingResult

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Staging results
# ---------------------------------------------------------------------------

def get_staging_result_by_experiment(
    experiment_id: int,
    db_path: str = "database.db",
) -> Optional[StagingResult]:
    """
    Fetch the most recent *valid* staging result for a given experiment.

    Returns
    -------
    StagingResult if found and is_valid=True, otherwise None.
    """
    query = """
    SELECT
        id,
        experiment_id,
        is_valid,
        match_type,
        identity,
        coverage,
        strand,
        frame,
        start_nt,
        end_nt_exclusive,
        wraps_origin,
        notes
    FROM staging_results
    WHERE experiment_id = ?
      AND is_valid = 1
    ORDER BY id DESC
    LIMIT 1;
    """

    with sqlite3.connect(db_path) as con:
        row = con.execute(query, (experiment_id,)).fetchone()

    if row is None:
        logger.warning(
            f"No valid staging result found for experiment_id={experiment_id}."
        )
        return None

    staging = StagingResult(
        id=row[0],
        experiment_id=row[1],
        is_valid=bool(row[2]),
        match_type=row[3],
        identity=row[4],
        coverage=row[5],
        strand=row[6],
        frame=row[7],
        start_nt=row[8],
        end_nt_exclusive=row[9],
        wraps_origin=bool(row[10]),
        notes=row[11],
    )

    logger.info(
        f"Staging result loaded: id={staging.id}, "
        f"strand={staging.strand}, frame={staging.frame}, "
        f"identity={staging.identity:.3f}, coverage={staging.coverage:.3f}"
    )

    return staging


# ---------------------------------------------------------------------------
# Variants
# ---------------------------------------------------------------------------

def get_variant_plasmid_sequence(
    plasmid_variant_index: str,
    db_path: str = "database.db",
) -> Optional[str]:
    """
    Fetch the assembled DNA sequence for a given variant from the variants table.
    Returns None if the variant is not found.
    """
    query = """
    SELECT gene_sequence
    FROM variants
    WHERE Plasmid_Variant_Index = ?
    LIMIT 1;
    """

    with sqlite3.connect(db_path) as con:
        row = con.execute(query, (plasmid_variant_index,)).fetchone()

    if row is None:
        logger.warning(f"No variant found for index={plasmid_variant_index}.")
        return None

    return row[0]


def get_wt_plasmid_sequence(
    experiment_id: int,
    db_path: str = "database.db",
) -> Optional[str]:
    """
    Fetch the WT plasmid sequence associated with an experiment.
    Stored in the experiments table as plasmid_sequence.
    Returns None if not found.
    """
    query = """
    SELECT plasmid_sequence
    FROM experiments
    WHERE id = ?
    LIMIT 1;
    """

    with sqlite3.connect(db_path) as con:
        row = con.execute(query, (experiment_id,)).fetchone()

    if row is None:
        logger.warning(f"No experiment found for id={experiment_id}.")
        return None

    return row[0]


def get_wt_protein_sequence(
    experiment_id: int,
    db_path: str = "database.db",
) -> Optional[str]:
    """
    Fetch the WT protein sequence for an experiment.
    Pulled from UniProt during staging and stored in uniprot_features.
    Returns None if not found.
    """
    query = """
    SELECT uf.features_json
    FROM uniprot_features uf
    WHERE uf.experiment_id = ?
    LIMIT 1;
    """

    with sqlite3.connect(db_path) as con:
        row = con.execute(query, (experiment_id,)).fetchone()

    if row is None:
        logger.warning(f"No UniProt features found for experiment_id={experiment_id}.")
        return None

    import json
    features = json.loads(row[0])
    return features.get("sequence", None)
