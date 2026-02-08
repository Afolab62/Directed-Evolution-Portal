import sqlite3
import numpy as np
import pandas as pd

abs_min = 1e-6 # absolute threshold
sc_min = 0.05 # 5% noise threshold

# the following function creates an SQLite table in the current working directory.

def database(db_path):
	'''
	Establishes a connection to a SQLite database file and creates a table
	to store activity score results. 
	
	Parameters
	----------
	db_path: str
		Pathway to the SQLite database file.
	
	Returns
	-------
	None

	'''
	
	# the following command, sqlite3.connect(), creates a connection to
	# the database in the working directory, implicitly creating one if 
	# it does not exist already. 
	
	con = sqlite3.connect('database.db')
	
	# the following command, con.cursor(), creates a cursor to execute
	# SQL statements and fetch results from SQL queries.
	
	cur = con.cursor()
	
	# the following command, cur.execute(...), executes the 'CREATE
	# TABLE' statement. 
        
	cur.execute('CREATE TABLE activity_scores(Directed_Evolution_Generation, DNA_Quantification_fg, Protein_Quantification_pg, dna_baseline, protein_baseline, dna_corrected, protein_corrected, activity_score, Control)')

	# the following command, con.commit(), commits transactions, 
	# a necessary step before changes can be saved in the database. 
	
	con.commit()
	con.close()

# the following function creates activity scores for each generation of a given variant.

def compute_activity_score(df):
	'''
	Computes a generation-normalised Activity Score.

	Activity Score is defined as baseline-corrected DNA yield divided by
	baseline-corrected protein yield:

	DNA* = DNA(g, v) - DNA(g, control)
	Protein* = Protein(g, v) - Protein(g, control)
	ActivityScore = DNA*/ Protein*

	Implementation Details:
	- DNA* is clipped at 0, so DNA is only counted above baseline.
	- A relative threshold is utilised, whereby there are two
	  minimal values and the larger is implemented:
	  The absolute minimal value (1e-6) prevents division by 
	  values close to zero. The other potential minimal value scales with 
	  the typical protein expression level -- 5% is taken as 
	  below this threshold is typically noise. 
          
	Parameters
	----------
 	df: DataFrame
		Input table containing DNA and Protein yield measurements
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
                        Baseline protein yield for the corresponding 
			generation.
                - dna_corrected: float
                        Baseline-corrected DNA yield (DNA*).
		- protein_corrected: float
			Baseline-corrected protein yield (Protein*).
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

	# the following command creates a copy of the dataframe. 

	out = df.copy()

	# the following commands convert the DNA and protein yield columns to
	# numeric values. this is to ensure the values are not registered as
	# string values.

	out['DNA_Quantification_fg'] = pd.to_numeric(out['DNA_Quantification_fg'], errors = 'coerce')
	out['Protein_Quantification_pg'] = pd.to_numeric(out['Protein_Quantification_pg'], errors = 'coerce')

	# the following commands verify the control experiment rows. a check-
	# point step ensures that there are control experiments present in the
	# data.

	controls = out[out['Control'] == True]
	if controls.empty:
		raise ValueError('No control data found')

	# the following commands create the baseline dna and protein yield values
	# by firstly grouping each generation's control experiments and then 
	# calculating the median value.	

	baselines = (
		controls.groupby('Directed_Evolution_Generation')
		.agg(
			dna_baseline = ('DNA_Quantification_fg', 'median'),
			protein_baseline = ('Protein_Quantification_pg', 'median'), 
		     )
	)

	# the following command attaches the baseline values to each row. 
	
	out = out.merge(baselines, on = 'Directed_Evolution_Generation', how = 'right')

	# the following commands calculate baseline-corrected DNA and protein 
	# yield. the implementation details are followed.

	out['dna_corrected'] = (out['DNA_Quantification_fg'] - out['dna_baseline']).clip(lower = 0)
	
	r_protein = out['Protein_Quantification_pg'] - out['protein_baseline']
	r_min = sc_min * out['protein_baseline']
	protein_min = np.maximum(abs_min, r_min)
	
	# the following command enforces the minimum allowed protein value.

	out['protein_corrected'] = np.where(r_protein < protein_min, protein_min, r_protein)
	
	# the following commands calculate the activity score.

	out['activity_score'] = out['dna_corrected']/ out['protein_corrected']
	out.loc[r_protein <= 0, 'activity_score'] = pd.NA

	return out

# the following command writes the activity scores to the SQLite database.

def write_activity_scores_to_db(scored_df, db_path):
	'''
	Writes computed activity score results to a SQLite database.

	Parameters
	----------
	scored_df: DataFrame	
		Dataframe containing activity score results produced by
		'compute_activity_score'.
	db_path: str
		Pathway to the SQLite database.

	Returns
	-------
	None	
	'''
	
	required_cols = {
                'Directed_Evolution_Generation',
                'DNA_Quantification_fg',
                'Protein_Quantification_pg',
                'Control',
		'dna_baseline',
		'protein_baseline',
		'dna_corrected',
		'protein_corrected',
		'activity_score'
        }

	# checkpoint

	missing = set(required_cols) - set(scored_df.columns)
	if missing:
		raise ValueError(f'Scored DataFrame Missing Required Columns : {sorted(missing)}')
	
	# the following command, .to_sql(),  opens a connection to the SQLite
	# database and writes the results table.
	
	with sqlite3.connect(db_path) as x:
		scored_df[required_cols].to_sql(
			'activity_scores',
			x,
			if_exists = if_exists,
			index = False
	)
