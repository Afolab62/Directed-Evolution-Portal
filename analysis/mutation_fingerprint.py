

import sqlite3
import pandas as pd


def reconstruct_lineage(selected_variant, db_path="database.db"):
    """
    Reconstruct the evolutionary lineage of a variant
    by following Parent_Plasmid_Variant backwards.

    Returns a DataFrame ordered chronologically
    (earliest generation first).
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

            variant_id = row[0]
            parent_id = row[1]
            generation = row[2]

            # Store current variant information
            lineage_records.append({
                "variant": variant_id,
                "parent": parent_id,
                "generation": generation
            })

            # Move to parent
            if parent_id is None:
                break

            current_variant = parent_id

    # Convert to DataFrame
    lineage_df = pd.DataFrame(lineage_records)

    # Sort chronologically
    lineage_df = lineage_df.sort_values(
        by="generation"
    ).reset_index(drop=True)

    return lineage_df

import sqlite3
import pandas as pd


def compute_mutation_fingerprint(
    selected_variant,
    db_path="database.db"
):
    """
    Reconstruct mutation fingerprint assuming mutations
    are stored per generation (non-cumulative model).

    Each variant contains only mutations introduced
    in that generation.
    """

    # Step 1: reconstruct lineage
    lineage_df = reconstruct_lineage(selected_variant, db_path)

    fingerprint_records = []

    mutation_query = """
    SELECT
        position,
        wt_aa,
        mut_aa
    FROM mutations
    WHERE Plasmid_Variant_Index = ?
    """

    with sqlite3.connect(db_path) as con:
        cursor = con.cursor()

        # Step 2: loop forward through lineage
        for _, row in lineage_df.iterrows():

            variant_id = row["variant"]
            generation = row["generation"]

            # Skip generation 0 (wildtype)
            if generation == 0:
                continue

            cursor.execute(mutation_query, (variant_id,))
            mutation_rows = cursor.fetchall()

            # Step 3: store mutations introduced in this generation
            for mut in mutation_rows:
                fingerprint_records.append({
                    "position": mut[0],
                    "wt_aa": mut[1],
                    "mut_aa": mut[2],
                    "generation_introduced": generation
                })

    fingerprint_df = pd.DataFrame(fingerprint_records)

    return fingerprint_df


import plotly.graph_objects as go
# fingerprint.py

import plotly.graph_objects as go


def plot_mutation_fingerprint(fingerprint_df, protein_length, variant_id):

    fig = go.Figure()

    # Protein backbone
    fig.add_shape(
        type="line",
        x0=1,
        x1=protein_length,
        y0=1,
        y1=1,
        line=dict(color="lightgrey", width=8)
    )

    # Plot mutations
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
                    f"Mutation: {row['wt_aa']}>{row['mut_aa']}<br>"
                    f"Generation: {row['generation_introduced']}"
                ),
                hoverinfo="text"
            )
        )

    fig.update_layout(
        title=f"Mutation Fingerprint — Variant {variant_id}",
        xaxis_title="Amino Acid Position",
        yaxis=dict(visible=False)
    )

    fig.show()








