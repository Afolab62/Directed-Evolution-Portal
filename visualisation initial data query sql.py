import sqlite3
import pandas as pd


def get_top_performers(db_path="database.db"):
    """
    Retrieve the top 10 performing variants from the final directed evolution generation.

    The function:
    - Excludes control variants
    - Restricts selection to the highest generation
    - Aggregates replicate activity scores (mean and replicate count)
    - Requires at least two replicates per variant
    - Computes standard deviation and standard error
    - Adds mutation summary statistics

    Returns
    -------
    pandas.DataFrame
        A dataframe containing variant-level performance and mutation summary.
    """

    with sqlite3.connect(db_path) as con:

       
        # 1. Variant ranking based on mean activity (replicate-level)
        # Aggregate replicate activity scores per variant.
        # Restrict to non-control variants from the final generation.
        # Enforce a minimum of two replicates for statistical robustness.
        ranking_query = """
        SELECT
            v.Plasmid_Variant_Index,
            v.Parent_Plasmid_Variant,
            v.Directed_Evolution_Generation,
            AVG(a.activity_score) AS Mean_Activity,
            COUNT(a.activity_score) AS Replicate_Count
        FROM variants v
        JOIN activity_scores a
            ON v.Plasmid_Variant_Index = a.Plasmid_Variant_Index
        WHERE v.Control = 0
          AND v.Directed_Evolution_Generation = (
                SELECT MAX(Directed_Evolution_Generation)
                FROM variants
          )
        GROUP BY v.Plasmid_Variant_Index
        HAVING COUNT(a.activity_score) >= 2
        ORDER BY Mean_Activity DESC
        LIMIT 10;
        """

        top_df = pd.read_sql(ranking_query, con)

        # If no qualifying variants are found, return empty dataframe
        if top_df.empty:
            return top_df

        # 2. Compute activity variability statistics (SD and SE)
        
        # Retrieve replicate-level activity data only for the selected
        # top-performing variants.
        top_variants = top_df["Plasmid_Variant_Index"].tolist()
        placeholders = ",".join(["?"] * len(top_variants))

        replicate_query = f"""
        SELECT
            Plasmid_Variant_Index,
            activity_score
        FROM activity_scores
        WHERE Plasmid_Variant_Index IN ({placeholders});
        """

        replicate_df = pd.read_sql(
            replicate_query,
            con,
            params=top_variants
        )

        # Calculate standard deviation of activity per variant
        sd_df = (
            replicate_df
            .groupby("Plasmid_Variant_Index")["activity_score"]
            .std()
            .reset_index()
            .rename(columns={"activity_score": "SD_Activity"})
        )

        # Merge SD into ranking dataframe
        top_df = top_df.merge(
            sd_df,
            on="Plasmid_Variant_Index",
            how="left"
        )

        # Calculate standard error (SD / sqrt(n))
        top_df["SE_Activity"] = (
            top_df["SD_Activity"] /
            top_df["Replicate_Count"] ** 0.5
        )

        
        # 3. Mutation summary statistics
       
        # Compute total, synonymous, and non-synonymous mutation counts
        # per variant independently of activity ranking.
        mutation_query = """
        SELECT
            Plasmid_Variant_Index,
            COUNT(*) AS Total_Mutations,
            SUM(CASE WHEN mutation_type = 'non-synonymous'
                     THEN 1 ELSE 0 END) AS Non_Synonymous_Count,
            SUM(CASE WHEN mutation_type = 'synonymous'
                     THEN 1 ELSE 0 END) AS Synonymous_Count
        FROM mutations
        GROUP BY Plasmid_Variant_Index;
        """

        mutation_df = pd.read_sql(mutation_query, con)

        # Merge mutation summary into ranking dataframe
        top_df = top_df.merge(
            mutation_df,
            on="Plasmid_Variant_Index",
            how="left"
        )

    # Ensure ranking order is preserved after merges
    top_df = (
        top_df
        .sort_values(by="Mean_Activity", ascending=False)
        .reset_index(drop=True)
    )

    return top_df











        

