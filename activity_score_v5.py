import sqlite3
import numpy as np
import pandas as pd

# activity scores for each generation of a given variant

def compute_activity_score(df):
	'''
	Computes a generation-normalised Activity Score.

	Activity Score is defined as baseline-corrected DNA yield divided by
	baseline-corrected protein concentration:

	DNA* = DNA(g, v) - DNA(g, control)
	Protein* = Protein(g, v) - Protein(g, control)
	ActivityScore = DNA*/ Protein*
	          
	Parameters
	----------
 	df: DataFrame
		Input table containing DNA and Protein concentration measurements
		for variants across generations. The dataframe must contain
		the following columns:
		- Plasmid_Variant_Index: str
			Unique Identifier for each variant.
		- Directed_Evolution_Generation: int
			Generation Identifier.
		- DNA_Quantification_fg: float
			Measured DNA yield for each variant.
		- Protein_Quantification_pg: float	
			Measured protein expression level yield for each 
			variant.
		- Control: bool
			True for baseline control measurements.
			False for variants.  

	Returns
	-------
	scored_df: DataFrame
		Copy of the input dataframe with additional columns:  
                - dna_baseline: float
                        Baseline DNA yield for the corresponding generation.
                - protein_baseline: float
                        Baseline protein concentration for the corresponding 
			generation.
                - dna*: float
                        Baseline-corrected DNA yield (DNA*).
		- protein*: float
			Baseline-corrected protein concentration (Protein*).
		- activity_score: float
			Generation-normalised Activity Score. 

	'''

	required = {
		'Directed_Evolution_Generation',
		'DNA_Quantification_fg', 
		'Protein_Quantification_pg', 		
		'Control'
	}
 	
	# checkpoint

	missing = required - set(df.columns)
	if missing:
		raise ValueError(f'Missing Required Columns: {sorted(missing)}')

	# copy of the dataframe 

	out = df.copy()

	# DNA yield and protein concentration columns converted to numeric values
	# (ensures the values are not registered as string values)

	out['DNA_Quantification_fg'] = pd.to_numeric(out['DNA_Quantification_fg'], errors = 'coerce')
	out['Protein_Quantification_pg'] = pd.to_numeric(out['Protein_Quantification_pg'], errors = 'coerce')

	# checkpoint

	controls = out[out['Control'] == True]
	if controls.empty:
		raise ValueError('No control data found')

	# baseline dna and protein yield values
	# (each generation's control experiments is grouped and then 
	# the median value calculated)	

	baselines = (
		controls.groupby('Directed_Evolution_Generation')
		.agg(
			dna_baseline = ('DNA_Quantification_fg', 'median'),
			protein_baseline = ('Protein_Quantification_pg', 'median'), 
		     )
	)

	# baseline values attached to each row 
	
	out = out.merge(baselines, on = 'Directed_Evolution_Generation', how = 'right')

	# baseline-corrected DNA and baseline-corrected protein yield calculation 

	out['dna*'] = (out['DNA_Quantification_fg'] - out['dna_baseline']).clip(lower = 0)
	
	out['protein*'] = out['Protein_Quantification_pg'] - out['protein_baseline']	
	
	# activity score calculation

	out['activity_score'] = out['dna*']/ out['protein*']
	out.loc[out['activity_score'] <= 0, 'activity_score'] = pd.NA

	return out

path = '/Users/minaal/Downloads/example.tsv'
df = pd.read_csv(path, sep = '\t')
x = compute_activity_score(df)
print(x)
                        
top_10 = x.sort_values(
        by = 'activity_score',  
        ascending = False
).head(10)      
        
print(top_10)   
        




