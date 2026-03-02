from __future__ import annotations

import logging
import sqlite3
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)


def load_variants_from_db(db_path: Path) -> dict[int, dict[str, Any]]:
    """
    Load the variant table from SQLite in the same shape used by
    mutation_analysis.load_variants_table().
    """
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    try:
        rows = conn.execute(
            """
            SELECT
                Plasmid_Variant_Index,
                Parent_Plasmid_Variant,
                Directed_Evolution_Generation,
                Assembled_DNA_Sequence,
                Protein_Sequence
            FROM activity_measurements
            ORDER BY Plasmid_Variant_Index
            """
        ).fetchall()
    finally:
        conn.close()

    variants: dict[int, dict[str, Any]] = {}
    for row in rows:
        try:
            variant_id = int(row["Plasmid_Variant_Index"])
        except (TypeError, ValueError):
            logger.warning("Skipping row with invalid Plasmid_Variant_Index: %r", row["Plasmid_Variant_Index"])
            continue

        parent = row["Parent_Plasmid_Variant"]
        generation = row["Directed_Evolution_Generation"]
        variants[variant_id] = {
            "variant_id": variant_id,
            "parent_variant_id": int(parent) if parent is not None else -1,
            "generation": int(generation) if generation is not None else -1,
            "dna_sequence": str(row["Assembled_DNA_Sequence"] or ""),
            "protein_sequence": str(row["Protein_Sequence"] or ""),
        }

    logger.info("Loaded %d variants from %s", len(variants), db_path)
    return variants


def get_all_variant_ids(db_path: Path) -> list[int]:
    """Return sorted variant ids from the activity_measurements table."""
    conn = sqlite3.connect(db_path)
    try:
        rows = conn.execute(
            """
            SELECT Plasmid_Variant_Index
            FROM activity_measurements
            WHERE Plasmid_Variant_Index IS NOT NULL
            ORDER BY Plasmid_Variant_Index
            """
        ).fetchall()
    finally:
        conn.close()

    return [int(row[0]) for row in rows]