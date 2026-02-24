"""
analysis/get_top_performers.py
--------------------------------
Query and return the top 10 performing variants by activity score.

This module only reads from the DB — all write operations happen
in the parsing/analysis pipeline, not here.
"""

import sqlite3
import pandas as pd
import logging

logger = logging.getLogger(__name__)


def get_top_performers(db_path: str = "database.db") -> pd.DataFrame:
    """
    Return the top 10 non-control variants ranked by activity_score (descending).

    Joins variants, activity_scores, and mutations tables.
    Excludes controls and rows with null activity scores.

    Parameters
    ----------
    db_path : str
        Path to the SQLite database.

    Returns
    -------
    pd.DataFrame with columns:
        Plasmid_Variant_Index, Parent_Plasmid_Variant,
        Directed_Evolution_Generation, DNA_Quantification_fg,
        Protein_Quantification_pg, activity_score,
        Total_Mutations, Synonymous_Count, Non_Synonymous_Count
    """
    query = """
    SELECT
        v.Plasmid_Variant_Index,
        v.Parent_Plasmid_Variant,
        v.Directed_Evolution_Generation,
        v.DNA_Quantification_fg,
        v.Protein_Quantification_pg,
        a.activity_score,
        m.Total_Mutations,
        m.Synonymous_Count,
        m.Non_Synonymous_Count

    FROM variants v

    JOIN activity_scores a
        ON v.Plasmid_Variant_Index = a.Plasmid_Variant_Index

    JOIN mutations m
        ON v.Plasmid_Variant_Index = m.Plasmid_Variant_Index

    WHERE v.Control = 0
      AND a.activity_score IS NOT NULL

    ORDER BY a.activity_score DESC

    LIMIT 10;
    """

    with sqlite3.connect(db_path) as con:
        df = pd.read_sql(query, con)

    if df.empty:
        logger.warning("No top performer data found — check that analysis has been run.")

    logger.info(f"Top performers fetched: {len(df)} rows.")
    return df
