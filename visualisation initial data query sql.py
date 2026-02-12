#This code connects to the SQLite database, executes your SQL query, loads the result into a Pandas DataFrame,
#  closes the connection automatically, and returns the DataFrame.
import sqlite3
import pandas as pd





def get_top_performers(db_path="database.db"):

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
        top_df = pd.read_sql(query, con)

    return top_df









        

