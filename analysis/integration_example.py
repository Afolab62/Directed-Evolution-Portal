"""
Integration Script: Sequence Analysis + Activity Metrics
Demonstrates how Section 4.a connects to Section 4.b

This script shows how to:
1. Load sequence analysis results from database
2. Load experimental quantification data
3. Combine them for activity score calculations
4. Prepare data for visualization (Section 4.c)

Author: BIO727P Group Project
Date: February 2026
"""

import sqlite3
import pandas as pd
import json
import numpy as np
from typing import Dict, List
import os


class IntegratedAnalysis:
    """
    Integration class combining sequence analysis with experimental metrics.
    Bridges Section 4.a (sequence) with Section 4.b (activity scoring).
    """
    
    def __init__(self, db_path: str, experimental_data_path: str):
        """
        Initialize with database and experimental data.
        
        Args:
            db_path: Path to sequence analysis database
            experimental_data_path: Path to TSV/JSON with quantification data
        """
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path)
        
        # Load experimental data
        if experimental_data_path.endswith('.tsv'):
            self.exp_data = pd.read_csv(experimental_data_path, sep='\t')
        elif experimental_data_path.endswith('.json'):
            self.exp_data = pd.read_json(experimental_data_path)
        else:
            raise ValueError("Experimental data must be TSV or JSON")
        
        # Load sequence analysis results
        self.variants_df = pd.read_sql_query(
            "SELECT * FROM variant_sequences", 
            self.conn
        )
        
        print(f"Loaded {len(self.variants_df)} variants from database")
        print(f"Loaded {len(self.exp_data)} experimental records")
    
    def merge_sequence_and_experimental_data(self) -> pd.DataFrame:
        """
        Merge sequence analysis with experimental quantification data.
        
        Returns:
            Combined DataFrame with all information
        """
        # Merge on variant index
        merged = self.variants_df.merge(
            self.exp_data,
            left_on='variant_index',
            right_on='Plasmid_Variant_Index',
            how='inner'
        )
        
        print(f"Merged dataset contains {len(merged)} records")
        return merged
    
    def calculate_activity_scores(self, merged_df: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate unified activity metric for each variant.
        
        Activity Score = (DNA_yield - baseline_DNA) / (Protein_yield - baseline_protein)
        
        This is a PLACEHOLDER implementation. Your team should define the exact
        activity score calculation based on project requirements.
        
        Args:
            merged_df: Merged sequence and experimental data
        
        Returns:
            DataFrame with added activity_score column
        """
        
        # Separate control samples (WT) from variants
        controls = merged_df[merged_df['Control'] == True].copy()
        variants = merged_df[merged_df['Control'] == False].copy()
        
        # Calculate baseline from controls (per generation)
        baseline_stats = controls.groupby('Directed_Evolution_Generation').agg({
            'DNA_Quantification_fg': ['mean', 'std'],
            'Protein_Quantification_pg': ['mean', 'std']
        })
        
        print("\nBaseline (WT Control) Statistics per Generation:")
        print(baseline_stats)
        
        # Calculate activity scores for variants
        # NOTE: This is an example calculation - adjust based on your methodology
        activity_scores = []
        
        for idx, row in variants.iterrows():
            gen = row['Directed_Evolution_Generation']
            
            # Get baseline for this generation
            if gen in baseline_stats.index:
                baseline_dna = baseline_stats.loc[gen, ('DNA_Quantification_fg', 'mean')]
                baseline_protein = baseline_stats.loc[gen, ('Protein_Quantification_pg', 'mean')]
            else:
                # Use overall mean if generation not in controls
                baseline_dna = controls['DNA_Quantification_fg'].mean()
                baseline_protein = controls['Protein_Quantification_pg'].mean()
            
            # Calculate normalized yields
            norm_dna = row['DNA_Quantification_fg'] / baseline_dna
            norm_protein = row['Protein_Quantification_pg'] / baseline_protein
            
            # Activity score: DNA yield normalized by protein expression
            # High DNA yield with low protein = efficient enzyme
            if norm_protein > 0:
                activity_score = norm_dna / norm_protein
            else:
                activity_score = 0
            
            activity_scores.append(activity_score)
        
        variants['activity_score'] = activity_scores
        
        # Combine back with controls (controls get baseline score of 1.0)
        controls['activity_score'] = 1.0
        
        result = pd.concat([variants, controls], ignore_index=True)
        result = result.sort_values(['Directed_Evolution_Generation', 'variant_index'])
        
        return result
    
    def identify_top_performers(self, df: pd.DataFrame, top_n: int = 10) -> pd.DataFrame:
        """
        Identify top N performing variants by activity score.
        
        Args:
            df: DataFrame with activity scores
            top_n: Number of top variants to return
        
        Returns:
            DataFrame with top performers
        """
        # Exclude controls from ranking
        variants_only = df[df['Control'] == False].copy()
        
        # Sort by activity score
        top_performers = variants_only.nlargest(top_n, 'activity_score')
        
        return top_performers
    
    def generate_activity_summary(self, df: pd.DataFrame) -> Dict:
        """
        Generate summary statistics about activity scores.
        
        Args:
            df: DataFrame with activity scores
        
        Returns:
            Dictionary with summary statistics
        """
        variants_only = df[df['Control'] == False].copy()
        
        summary = {
            'total_variants': len(variants_only),
            'mean_activity': variants_only['activity_score'].mean(),
            'median_activity': variants_only['activity_score'].median(),
            'std_activity': variants_only['activity_score'].std(),
            'max_activity': variants_only['activity_score'].max(),
            'min_activity': variants_only['activity_score'].min(),
            'num_above_baseline': (variants_only['activity_score'] > 1.0).sum(),
            'num_below_baseline': (variants_only['activity_score'] < 1.0).sum(),
        }
        
        return summary
    
    def prepare_visualization_data(self, df: pd.DataFrame) -> Dict:
        """
        Prepare data structures for Section 4.c visualizations.
        
        Args:
            df: Complete DataFrame with all metrics
        
        Returns:
            Dictionary with data formatted for different plots
        """
        
        viz_data = {}
        
        # 1. Activity distribution per generation (for box plots)
        viz_data['activity_by_generation'] = df.groupby('Directed_Evolution_Generation')['activity_score'].apply(list).to_dict()
        
        # 2. Top performers table data
        top_10 = self.identify_top_performers(df, top_n=10)
        viz_data['top_performers'] = top_10[[
            'variant_index',
            'generation',
            'activity_score',
            'num_mutations',
            'num_synonymous',
            'num_nonsynonymous',
            'DNA_Quantification_fg',
            'Protein_Quantification_pg'
        ]].to_dict('records')
        
        # 3. Mutation vs activity correlation
        viz_data['mutation_activity_data'] = df[[
            'variant_index',
            'num_mutations',
            'num_nonsynonymous',
            'activity_score',
            'generation'
        ]].to_dict('records')
        
        # 4. Generation progression
        gen_summary = df.groupby('Directed_Evolution_Generation').agg({
            'activity_score': ['mean', 'std', 'max'],
            'num_mutations': 'mean'
        })
        gen_summary.columns = ['_'.join(col).strip('_') for col in gen_summary.columns.values]
        gen_summary = gen_summary.reset_index()
        viz_data['generation_progression'] = gen_summary.to_dict('records')
        
        return viz_data
    
    def export_integrated_results(self, output_dir: str):
        """
        Export all integrated analysis results.
        
        Args:
            output_dir: Directory for output files
        """
        os.makedirs(output_dir, exist_ok=True)
        
        # Merge and calculate
        merged = self.merge_sequence_and_experimental_data()
        with_scores = self.calculate_activity_scores(merged)
        
        # Export complete dataset
        output_file = os.path.join(output_dir, "integrated_analysis.csv")
        with_scores.to_csv(output_file, index=False)
        print(f"\nExported integrated data to: {output_file}")
        
        # Export top performers
        top_performers = self.identify_top_performers(with_scores, top_n=10)
        top_file = os.path.join(output_dir, "top_10_performers.csv")
        top_performers.to_csv(top_file, index=False)
        print(f"Exported top performers to: {top_file}")
        
        # Export visualization-ready data
        viz_data = self.prepare_visualization_data(with_scores)
        viz_file = os.path.join(output_dir, "visualization_data.json")
        with open(viz_file, 'w') as f:
            json.dump(viz_data, f, indent=2)
        print(f"Exported visualization data to: {viz_file}")
        
        # Print summary
        summary = self.generate_activity_summary(with_scores)
        print("\n" + "="*80)
        print("ACTIVITY ANALYSIS SUMMARY")
        print("="*80)
        for key, value in summary.items():
            if isinstance(value, float):
                print(f"{key}: {value:.4f}")
            else:
                print(f"{key}: {value}")
        
        print("\nTop 10 Performers:")
        print("-"*80)
        print(top_performers[['variant_index', 'generation', 'activity_score', 'num_mutations']].to_string(index=False))
        
        return with_scores
    
    def close(self):
        """Close database connection"""
        if self.conn:
            self.conn.close()


def main():
    """
    Main function demonstrating integrated analysis workflow.
    """
    
    print("="*80)
    print("INTEGRATED ANALYSIS: Sequence + Activity Metrics")
    print("Connecting Section 4.a → 4.b → 4.c")
    print("="*80)
    
    # File paths
    DB_PATH = "/mnt/user-data/outputs/directed_evolution_analysis.db"
    EXP_DATA_PATH = "/mnt/user-data/uploads/DE_BSU_Pol_Batch_1.tsv"
    OUTPUT_DIR = "/mnt/user-data/outputs/integrated_analysis"
    
    # Check if database exists
    if not os.path.exists(DB_PATH):
        print(f"\nERROR: Database not found at {DB_PATH}")
        print("Please run process_data.py first!")
        return
    
    # Initialize integrated analysis
    print(f"\nInitializing integrated analysis...")
    analysis = IntegratedAnalysis(DB_PATH, EXP_DATA_PATH)
    
    # Run complete analysis and export
    print("\nRunning integrated analysis and exporting results...")
    integrated_df = analysis.export_integrated_results(OUTPUT_DIR)
    
    # Example: Access specific variant details
    print("\n" + "="*80)
    print("EXAMPLE: Detailed Analysis of a Specific Variant")
    print("="*80)
    
    # Get variant with highest activity score
    best_variant = integrated_df[integrated_df['Control'] == False].nlargest(1, 'activity_score').iloc[0]
    
    print(f"\nBest Performing Variant: {int(best_variant['variant_index'])}")
    print(f"  Generation: {int(best_variant['generation'])}")
    print(f"  Activity Score: {best_variant['activity_score']:.4f}")
    print(f"  Total Mutations: {int(best_variant['num_mutations'])}")
    print(f"  Non-synonymous: {int(best_variant['num_nonsynonymous'])}")
    print(f"  DNA Yield: {best_variant['DNA_Quantification_fg']:.2f} fg")
    print(f"  Protein Yield: {best_variant['Protein_Quantification_pg']:.2f} pg")
    
    # Close connection
    analysis.close()
    
    print("\n" + "="*80)
    print("INTEGRATION COMPLETE")
    print("="*80)
    print(f"\nAll results saved to: {OUTPUT_DIR}")
    print("\nNext Steps:")
    print("1. Use 'integrated_analysis.csv' for further analysis")
    print("2. Use 'visualization_data.json' for Section 4.c plots")
    print("3. Use 'top_10_performers.csv' for detailed variant characterization")


if __name__ == "__main__":
    main()
