import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import sqlite3


def get_activity_scores_by_generation(db_path="database.db") -> pd.DataFrame:
    """
    Query the database for activity scores grouped by generation.
    Excludes controls and null activity scores.

    Returns a DataFrame with columns:
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

    return df


def plot_activity_score_distribution(db_path="database.db") -> None:
    """
    Plot the distribution of activity scores for each generation
    as a violin + strip plot.
    """
    df = get_activity_scores_by_generation(db_path)

    if df.empty:
        print("No activity score data found.")
        return

    fig, ax = plt.subplots(figsize=(10, 6))

    sns.violinplot(
        data=df,
        x="Directed_Evolution_Generation",
        y="activity_score",
        inner=None,
        color="lightblue",
        ax=ax
    )

    sns.stripplot(
        data=df,
        x="Directed_Evolution_Generation",
        y="activity_score",
        color="black",
        alpha=0.5,
        size=4,
        ax=ax
    )

    ax.set_title("Activity Score Distribution by Generation")
    ax.set_xlabel("Directed Evolution Generation")
    ax.set_ylabel("Activity Score")

    plt.tight_layout()
    plt.show()
