'''
Date: 2026-01-27
Here I perform an initial exloration of sample data 
1)The TSV file hads in total 10 rows and 7 columns 
Columns name:Plasmid_Variant_Index', 'Parent_Plasmid_Variant',
       'Directed_Evolution_Generation', 'Assembled_DNA_Sequence',
       'DNA_Quantification_fg', 'Protein_Quantification_pg', 'Control'],
Numeric columns(['Plasmid_Variant_Index', 'Parent_Plasmid_Variant',
       'Directed_Evolution_Generation', 'DNA_Quantification_fg',
       'Protein_Quantification_pg'],

Text columns(['Assembled_DNA_Sequence']
boolean=Control column


2)CSV
3)fasta.file
'''
#1
import pandas as pd
df_tsv=pd.read_csv("/Users/virginialaspina/Desktop/group_project/Directed-Evolution-Portal/data/sequences/DE_BSU_Pol_Batch_1.tsv" ,delimiter='\t')
print(df_tsv.head(10))
print(len(df_tsv))
print(df_tsv.columns)
numeric_cols = df_tsv.select_dtypes(include="number").columns
text_cols =df_tsv.select_dtypes(include ="object").columns
print(f"Numeric columns{numeric_cols}")
print(f"Text columns{text_cols}")
print(df_tsv.select_dtypes(include="boolean").columns)



#2

