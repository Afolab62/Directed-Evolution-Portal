"""
analysis/activity_score_distribution_plot.py
---------------------------------------------
Generate a violin + strip plot of activity score distributions
grouped by directed evolution generation.

Used in the web portal to let users quickly judge whether
activity is increasing across generations.
"""

import sqlite3
import logging

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

logger = logging.getLogger(__name__)


def get_activity_scores_by_generation(db_path: str = "database.db") -> pd.DataFrame:
    """
    Query activity scores grouped by generation.
    Excludes control samples and rows with null activity scores.

    Parameters
    ----------
    db_path : str
        Path to the SQLite database.

    Returns
    -------
    pd.DataFrame with columns:
        Directed_Evolution_Generation, activity_score
    """
    query = """
    SELECT
        v.Directed_Evolution_Generation,
        a.activity_score
    FROM variants v
    JOIN activity_scores a
        ON v.Plasmid_Variant_Index = a.Plasmid_Variant_Index
    WHERE v.Control = 0
      AND a.activity_score IS NOT NULL
    ORDER BY v.Directed_Evolution_Generation;
    """

    with sqlite3.connect(db_path) as con:
        df = pd.read_sql_query(query, con)

    logger.info(
        f"Activity scores fetched: {len(df)} rows across "
        f"{df['Directed_Evolution_Generation'].nunique()} generations."
    )
    return df


def plot_activity_score_distribution(
    db_path: str = "database.db",
    save_path: str = None,
) -> plt.Figure:
    """
    Plot the distribution of activity scores for each generation
    as a violin + strip plot.

    Parameters
    ----------
    db_path : str
        Path to the SQLite database.
    save_path : str, optional
        If provided, save the figure to this path instead of displaying it.
        Useful when generating PDFs or serving images via Flask.

    Returns
    -------
    matplotlib Figure object (for further customisation or embedding).
    """
    df = get_activity_scores_by_generation(db_path)

    if df.empty:
        logger.warning("No activity score data found — cannot generate plot.")
        return None

    fig, ax = plt.subplots(figsize=(10, 6))

    sns.violinplot(
        data=df,
        x="Directed_Evolution_Generation",
        y="activity_score",
        inner=None,
        color="lightblue",
        ax=ax,
    )

    sns.stripplot(
        data=df,
        x="Directed_Evolution_Generation",
        y="activity_score",
        color="black",
        alpha=0.5,
        size=4,
        ax=ax,
    )

    ax.set_title("Activity Score Distribution by Generation")
    ax.set_xlabel("Directed Evolution Generation")
    ax.set_ylabel("Activity Score")

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        logger.info(f"Activity score distribution plot saved to {save_path}.")
    else:
        plt.show()

    return fig
