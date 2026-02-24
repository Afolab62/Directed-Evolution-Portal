"""
analysis/mutation_fingerprint.py
---------------------------------
Reconstruct the evolutionary lineage of a variant and build its
mutation fingerprint — the cumulative set of amino acid changes
introduced across all generations leading to that variant.

Used for the optional "Mutation Fingerprinting by Generation" visualisation
described in the project brief.
"""

import sqlite3
import logging
from typing import Optional

import pandas as pd
import plotly.graph_objects as go

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Lineage reconstruction
# ---------------------------------------------------------------------------

def reconstruct_lineage(
    selected_variant: str,
    db_path: str = "database.db",
) -> pd.DataFrame:
    """
    Reconstruct the evolutionary lineage of a variant by following
    Parent_Plasmid_Variant backwards through the variants table.

    Parameters
    ----------
    selected_variant : str
        Plasmid_Variant_Index of the variant to trace.
    db_path : str
        Path to the SQLite database.

    Returns
    -------
    pd.DataFrame ordered chronologically (earliest generation first),
    with columns: variant, parent, generation.
    """
    lineage_records = []
    current_variant = selected_variant

    query = """
    SELECT
        Plasmid_Variant_Index,
        Parent_Plasmid_Variant,
        Directed_Evolution_Generation
    FROM variants
    WHERE Plasmid_Variant_Index = ?
    """

    with sqlite3.connect(db_path) as con:
        cursor = con.cursor()

        while current_variant is not None:
            cursor.execute(query, (current_variant,))
            row = cursor.fetchone()

            if row is None:
                break

            lineage_records.append({
                "variant": row[0],
                "parent": row[1],
                "generation": row[2],
            })

            # Walk backwards to parent; stop at root (parent is NULL)
            current_variant = row[1]

    if not lineage_records:
        logger.warning(f"No lineage found for variant={selected_variant}.")
        return pd.DataFrame(columns=["variant", "parent", "generation"])

    lineage_df = (
        pd.DataFrame(lineage_records)
        .sort_values("generation")
        .reset_index(drop=True)
    )

    logger.info(
        f"Lineage for {selected_variant}: "
        f"{len(lineage_df)} ancestors across "
        f"{lineage_df['generation'].max()} generation(s)."
    )
    return lineage_df


# ---------------------------------------------------------------------------
# Mutation fingerprint computation
# ---------------------------------------------------------------------------

def compute_mutation_fingerprint(
    selected_variant: str,
    db_path: str = "database.db",
) -> pd.DataFrame:
    """
    Compute the cumulative mutation fingerprint for a variant.

    Assumes mutations in the DB are stored per-generation (non-cumulative):
    each variants row contains only the mutations introduced *in that generation*.
    The fingerprint is built by walking forward through the lineage and collecting
    all mutations introduced at each step.

    Parameters
    ----------
    selected_variant : str
        Plasmid_Variant_Index of the variant to fingerprint.
    db_path : str
        Path to the SQLite database.

    Returns
    -------
    pd.DataFrame with columns:
        position, wt_aa, mut_aa, generation_introduced
    Empty DataFrame if no mutations found.
    """
    lineage_df = reconstruct_lineage(selected_variant, db_path)

    if lineage_df.empty:
        return pd.DataFrame(columns=["position", "wt_aa", "mut_aa", "generation_introduced"])

    mutation_query = """
    SELECT
        position,
        wt_aa,
        mut_aa
    FROM mutations
    WHERE Plasmid_Variant_Index = ?
    """

    fingerprint_records = []

    with sqlite3.connect(db_path) as con:
        cursor = con.cursor()

        for _, row in lineage_df.iterrows():
            variant_id = row["variant"]
            generation = row["generation"]

            # Skip generation 0 (wildtype — no mutations by definition)
            if generation == 0:
                continue

            cursor.execute(mutation_query, (variant_id,))
            for position, wt_aa, mut_aa in cursor.fetchall():
                fingerprint_records.append({
                    "position": position,
                    "wt_aa": wt_aa,
                    "mut_aa": mut_aa,
                    "generation_introduced": generation,
                })

    fingerprint_df = pd.DataFrame(fingerprint_records)

    logger.info(
        f"Mutation fingerprint for {selected_variant}: "
        f"{len(fingerprint_df)} total mutations across lineage."
    )
    return fingerprint_df


# ---------------------------------------------------------------------------
# Visualisation
# ---------------------------------------------------------------------------

def plot_mutation_fingerprint(
    fingerprint_df: pd.DataFrame,
    protein_length: int,
    variant_id: str,
    save_path: Optional[str] = None,
) -> go.Figure:
    """
    Plot the mutation fingerprint as a linear protein backbone with
    mutation markers coloured by generation introduced.

    Parameters
    ----------
    fingerprint_df : pd.DataFrame
        Output of compute_mutation_fingerprint().
    protein_length : int
        Total length of the WT protein (number of amino acids).
    variant_id : str
        Label for the plot title.
    save_path : str, optional
        If provided, save the figure as HTML to this path.

    Returns
    -------
    plotly Figure object.
    """
    fig = go.Figure()

    # Draw protein backbone as a grey line
    fig.add_shape(
        type="line",
        x0=1,
        x1=protein_length,
        y0=1,
        y1=1,
        line=dict(color="lightgrey", width=8),
    )

    if fingerprint_df.empty:
        logger.warning(f"No mutations to plot for variant={variant_id}.")
    else:
        for _, row in fingerprint_df.iterrows():
            fig.add_trace(
                go.Scatter(
                    x=[row["position"]],
                    y=[1.05],
                    mode="markers",
                    marker=dict(size=12),
                    name=f"Gen {row['generation_introduced']}",
                    hovertext=(
                        f"Position: {row['position']}<br>"
                        f"Mutation: {row['wt_aa']}→{row['mut_aa']}<br>"
                        f"Generation introduced: {row['generation_introduced']}"
                    ),
                    hoverinfo="text",
                )
            )

    fig.update_layout(
        title=f"Mutation Fingerprint — Variant {variant_id}",
        xaxis_title="Amino Acid Position",
        xaxis=dict(range=[0, protein_length + 10]),
        yaxis=dict(visible=False),
        legend_title="Generation",
    )

    if save_path:
        fig.write_html(save_path)
        logger.info(f"Mutation fingerprint saved to {save_path}.")
    else:
        fig.show()

    return fig
