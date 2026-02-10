"""
Database Query Utility for Directed Evolution Analysis
Provides functions to query and export analysis results

Author: BIO727P Group Project
Date: February 2026
"""

import sqlite3
import csv
import json
from typing import List, Dict, Optional
import pandas as pd


class AnalysisQueryTool:
    """Tool for querying and exporting analysis results from the database"""
    
    def __init__(self, db_path: str = "directed_evolution_analysis.db"):
        """
        Initialize query tool.
        
        Args:
            db_path: Path to SQLite database
        """
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path)
        self.conn.row_factory = sqlite3.Row  # Enable column access by name
    
    def get_variant_details(self, variant_index: int) -> Optional[Dict]:
        """
        Get complete details for a specific variant.
        
        Args:
            variant_index: Variant identifier
        
        Returns:
            Dictionary with variant details including mutations
        """
        cursor = self.conn.cursor()
        
        # Get variant info
        cursor.execute("""
            SELECT * FROM variant_sequences WHERE variant_index = ?
        """, (variant_index,))
        
        row = cursor.fetchone()
        if not row:
            return None
        
        variant = dict(row)
        
        # Get mutations
        cursor.execute("""
            SELECT * FROM mutations 
            WHERE variant_index = ? 
            ORDER BY position
        """, (variant_index,))
        
        mutations = [dict(row) for row in cursor.fetchall()]
        variant['mutations'] = mutations
        
        return variant
    
    def get_top_variants_by_generation(self, generation: int, top_n: int = 10) -> List[Dict]:
        """
        Get top N variants from a specific generation.
        Note: This returns variants sorted by variant_index. 
        For actual top performers, you need activity scores from Section 4.b
        
        Args:
            generation: Generation number
            top_n: Number of top variants to return
        
        Returns:
            List of variant dictionaries
        """
        cursor = self.conn.cursor()
        
        cursor.execute("""
            SELECT * FROM variant_sequences 
            WHERE generation = ?
            ORDER BY num_mutations ASC
            LIMIT ?
        """, (generation, top_n))
        
        return [dict(row) for row in cursor.fetchall()]
    
    def get_mutation_frequency_by_position(self) -> pd.DataFrame:
        """
        Calculate mutation frequency at each amino acid position.
        
        Returns:
            DataFrame with position and mutation counts
        """
        query = """
            SELECT 
                position,
                wt_aa,
                mut_aa,
                COUNT(*) as frequency,
                mutation_type
            FROM mutations
            GROUP BY position, wt_aa, mut_aa, mutation_type
            ORDER BY frequency DESC, position
        """
        
        return pd.read_sql_query(query, self.conn)
    
    def get_generation_summary(self) -> pd.DataFrame:
        """
        Get summary statistics for each generation.
        
        Returns:
            DataFrame with generation-level statistics
        """
        query = """
            SELECT 
                generation,
                COUNT(*) as num_variants,
                AVG(num_mutations) as avg_mutations,
                AVG(num_synonymous) as avg_synonymous,
                AVG(num_nonsynonymous) as avg_nonsynonymous,
                MAX(num_mutations) as max_mutations,
                MIN(num_mutations) as min_mutations
            FROM variant_sequences
            GROUP BY generation
            ORDER BY generation
        """
        
        return pd.read_sql_query(query, self.conn)
    
    def get_all_variants(self) -> pd.DataFrame:
        """
        Get all variant sequences.
        
        Returns:
            DataFrame with all variants
        """
        query = """
            SELECT * FROM variant_sequences
            ORDER BY generation, variant_index
        """
        
        return pd.read_sql_query(query, self.conn)
    
    def get_all_mutations(self) -> pd.DataFrame:
        """
        Get all mutations across all variants.
        
        Returns:
            DataFrame with all mutations
        """
        query = """
            SELECT m.*, v.generation
            FROM mutations m
            JOIN variant_sequences v ON m.variant_index = v.variant_index
            ORDER BY v.generation, m.variant_index, m.position
        """
        
        return pd.read_sql_query(query, self.conn)
    
    def export_variants_to_csv(self, output_path: str):
        """
        Export all variant data to CSV file.
        
        Args:
            output_path: Path for output CSV file
        """
        df = self.get_all_variants()
        df.to_csv(output_path, index=False)
        print(f"Exported {len(df)} variants to {output_path}")
    
    def export_mutations_to_csv(self, output_path: str):
        """
        Export all mutation data to CSV file.
        
        Args:
            output_path: Path for output CSV file
        """
        df = self.get_all_mutations()
        df.to_csv(output_path, index=False)
        print(f"Exported {len(df)} mutations to {output_path}")
    
    def export_variant_with_mutations_to_json(self, variant_index: int, output_path: str):
        """
        Export a single variant with all its mutations to JSON.
        
        Args:
            variant_index: Variant to export
            output_path: Path for output JSON file
        """
        variant = self.get_variant_details(variant_index)
        
        if variant:
            with open(output_path, 'w') as f:
                json.dump(variant, f, indent=2)
            print(f"Exported variant {variant_index} to {output_path}")
        else:
            print(f"Variant {variant_index} not found in database")
    
    def get_lineage_trace(self, variant_index: int) -> List[Dict]:
        """
        Trace the lineage of a variant back to the WT (generation 0).
        
        Args:
            variant_index: Variant to trace
        
        Returns:
            List of variant dictionaries in lineage order (WT to current)
        """
        lineage = []
        current_index = variant_index
        
        while current_index is not None and current_index >= 0:
            cursor = self.conn.cursor()
            cursor.execute("""
                SELECT * FROM variant_sequences WHERE variant_index = ?
            """, (current_index,))
            
            row = cursor.fetchone()
            if not row:
                break
            
            variant = dict(row)
            lineage.insert(0, variant)  # Add to beginning
            
            current_index = variant['parent_variant']
        
        return lineage
    
    def get_mutation_accumulation(self, variant_index: int) -> Dict:
        """
        Get cumulative mutations along a variant's lineage.
        
        Args:
            variant_index: Variant to analyze
        
        Returns:
            Dictionary with mutation accumulation data
        """
        lineage = self.get_lineage_trace(variant_index)
        
        accumulation = {
            'lineage_length': len(lineage),
            'generations': [],
            'cumulative_mutations': [],
            'mutations_per_generation': []
        }
        
        for variant in lineage:
            accumulation['generations'].append(variant['generation'])
            accumulation['cumulative_mutations'].append(variant['num_mutations'])
            
            # Calculate new mutations in this generation
            if len(accumulation['mutations_per_generation']) == 0:
                new_muts = variant['num_mutations']
            else:
                new_muts = variant['num_mutations'] - accumulation['cumulative_mutations'][-2]
            
            accumulation['mutations_per_generation'].append(new_muts)
        
        return accumulation
    
    def close(self):
        """Close database connection"""
        if self.conn:
            self.conn.close()


def generate_summary_report(db_path: str, output_dir: str = "/mnt/user-data/outputs"):
    """
    Generate comprehensive summary report with multiple exports.
    
    Args:
        db_path: Path to database
        output_dir: Directory for output files
    """
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    print("=" * 80)
    print("Generating Analysis Summary Report")
    print("=" * 80)
    
    tool = AnalysisQueryTool(db_path)
    
    # 1. Export all variants
    print("\n1. Exporting all variants...")
    variants_csv = os.path.join(output_dir, "all_variants.csv")
    tool.export_variants_to_csv(variants_csv)
    
    # 2. Export all mutations
    print("\n2. Exporting all mutations...")
    mutations_csv = os.path.join(output_dir, "all_mutations.csv")
    tool.export_mutations_to_csv(mutations_csv)
    
    # 3. Generation summary
    print("\n3. Generating generation summary...")
    gen_summary = tool.get_generation_summary()
    gen_summary_csv = os.path.join(output_dir, "generation_summary.csv")
    gen_summary.to_csv(gen_summary_csv, index=False)
    print(f"   Saved to {gen_summary_csv}")
    
    # 4. Mutation frequency analysis
    print("\n4. Analyzing mutation frequencies...")
    mut_freq = tool.get_mutation_frequency_by_position()
    mut_freq_csv = os.path.join(output_dir, "mutation_frequencies.csv")
    mut_freq.to_csv(mut_freq_csv, index=False)
    print(f"   Saved to {mut_freq_csv}")
    
    # 5. Print summary statistics
    print("\n5. Summary Statistics:")
    print("-" * 80)
    all_variants = tool.get_all_variants()
    print(f"Total variants analyzed: {len(all_variants)}")
    print(f"Generations: {all_variants['generation'].min()} to {all_variants['generation'].max()}")
    print(f"Average mutations per variant: {all_variants['num_mutations'].mean():.2f}")
    print(f"Total synonymous mutations: {all_variants['num_synonymous'].sum()}")
    print(f"Total non-synonymous mutations: {all_variants['num_nonsynonymous'].sum()}")
    
    print("\nGeneration Breakdown:")
    print(gen_summary.to_string(index=False))
    
    tool.close()
    
    print("\n" + "=" * 80)
    print(f"All reports saved to: {output_dir}")
    print("=" * 80)


def main():
    """Example usage of the query tool"""
    
    db_path = "/mnt/user-data/outputs/directed_evolution_analysis.db"
    
    # Check if database exists
    import os
    if not os.path.exists(db_path):
        print(f"Database not found at {db_path}")
        print("Please run process_data.py first to generate the database.")
        return
    
    # Generate comprehensive report
    generate_summary_report(db_path)
    
    # Example: Query specific variant
    print("\n\nExample: Querying variant details...")
    tool = AnalysisQueryTool(db_path)
    
    variant = tool.get_variant_details(1)
    if variant:
        print(f"\nVariant {variant['variant_index']}:")
        print(f"  Generation: {variant['generation']}")
        print(f"  Parent: {variant['parent_variant']}")
        print(f"  Total mutations: {variant['num_mutations']}")
        print(f"  Synonymous: {variant['num_synonymous']}")
        print(f"  Non-synonymous: {variant['num_nonsynonymous']}")
        print(f"  Protein sequence length: {len(variant['protein_sequence'])}")
        
        if variant['mutations']:
            print(f"\n  First 5 mutations:")
            for mut in variant['mutations'][:5]:
                print(f"    {mut['wt_aa']}{mut['position']}{mut['mut_aa']} ({mut['mutation_type']})")
    
    tool.close()


if __name__ == "__main__":
    main()
