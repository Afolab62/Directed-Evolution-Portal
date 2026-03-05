"""
Service for calculating activity scores for experimental variants with corrected edge case handling.

Activity Score = (DNA yield / DNA baseline) / (Protein yield / Protein baseline)

This fold-change ratio approach is more robust than simple subtraction and
naturally handles edge cases (negative values, division by zero).
"""

import pandas as pd
import numpy as np
from typing import Dict, List


class ActivityScoreCalculator:
    """Calculator for generation-normalized activity scores"""
    
    def __init__(self, epsilon: float = 0.01):
        """
        Initialize calculator
        
        Args:
            epsilon: Small value to prevent division by zero (default: 0.01)
        """
        self.epsilon = epsilon
    
    def calculate_baselines(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate baseline DNA and protein yields per generation from controls
        
        Args:
            df: DataFrame with columns: generation, dna_yield, protein_yield, is_control
        
        Returns:
            DataFrame with columns: generation, dna_baseline, protein_baseline
        
        Raises:
            ValueError: If no control data is found
        """
        # Filter controls
        controls = df[df['is_control'] == True].copy()
        
        if controls.empty:
            raise ValueError("No control data found. Controls are required for activity score calculation.")
        
        # Check for any generation without controls
        all_generations = set(df['generation'].unique())
        control_generations = set(controls['generation'].unique())
        missing_controls = all_generations - control_generations
        
        if missing_controls:
            print(f"Warning: Generations {sorted(missing_controls)} have no control data. "
                  f"Using overall median as baseline for these generations.")
        
        # Calculate baselines per generation (median is robust to outliers)
        baselines = controls.groupby('generation').agg(
            dna_baseline=('dna_yield', 'median'),
            protein_baseline=('protein_yield', 'median')
        ).reset_index()
        
        # For generations without controls, use overall median
        if missing_controls:
            overall_dna_baseline = controls['dna_yield'].median()
            overall_protein_baseline = controls['protein_yield'].median()
            
            for gen in missing_controls:
                baselines = pd.concat([
                    baselines,
                    pd.DataFrame({
                        'generation': [gen],
                        'dna_baseline': [overall_dna_baseline],
                        'protein_baseline': [overall_protein_baseline]
                    })
                ], ignore_index=True)
        
        return baselines
    
    def calculate_activity_scores(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate activity scores for all variants
        
        Activity Score = (DNA / DNA_baseline) / (Protein / Protein_baseline)
        
        This represents the fold-change in DNA yield normalized by fold-change
        in protein expression. High scores = efficient enzymes (more DNA per unit protein).
        
        Args:
            df: DataFrame with columns: generation, dna_yield, protein_yield, is_control
        
        Returns:
            DataFrame with added columns: dna_baseline, protein_baseline, activity_score
        
        Raises:
            ValueError: If required columns are missing or no controls found
        """
        # Validate required columns
        required_cols = {'generation', 'dna_yield', 'protein_yield', 'is_control'}
        missing_cols = required_cols - set(df.columns)
        if missing_cols:
            raise ValueError(f"Missing required columns: {sorted(missing_cols)}")
        
        # Create output dataframe
        result = df.copy()
        
        # Ensure numeric types
        result['dna_yield'] = pd.to_numeric(result['dna_yield'], errors='coerce')
        result['protein_yield'] = pd.to_numeric(result['protein_yield'], errors='coerce')
        result['generation'] = pd.to_numeric(result['generation'], errors='coerce').astype('Int64')
        
        # Calculate baselines
        baselines = self.calculate_baselines(result)
        
        # Merge baselines with data
        result = result.merge(baselines, on='generation', how='left')
        
        # Apply minimum threshold to prevent division by zero/near-zero
        result['dna_baseline_safe'] = result['dna_baseline'].clip(lower=self.epsilon)
        result['protein_baseline_safe'] = result['protein_baseline'].clip(lower=self.epsilon)
        result['dna_yield_safe'] = result['dna_yield'].clip(lower=self.epsilon)
        result['protein_yield_safe'] = result['protein_yield'].clip(lower=self.epsilon)
        
        # Calculate fold-change ratios for non-control samples
        non_control_mask = result['is_control'] == False
        
        result['activity_score'] = np.nan
        result.loc[non_control_mask, 'activity_score'] = (
            (result.loc[non_control_mask, 'dna_yield_safe'] / 
             result.loc[non_control_mask, 'dna_baseline_safe']) /
            (result.loc[non_control_mask, 'protein_yield_safe'] / 
             result.loc[non_control_mask, 'protein_baseline_safe'])
        )
        
        # Controls should have NA activity scores (not meaningful)
        result.loc[~non_control_mask, 'activity_score'] = np.nan
        
        # Clean up temporary columns
        result = result.drop(columns=[
            'dna_baseline_safe', 
            'protein_baseline_safe',
            'dna_yield_safe',
            'protein_yield_safe'
        ])
        
        return result
    
    def get_top_performers(
        self, 
        df: pd.DataFrame, 
        n: int = 10,
        exclude_controls: bool = True
    ) -> pd.DataFrame:
        """
        Get top N variants by activity score
        
        Args:
            df: DataFrame with activity_score column
            n: Number of top performers to return
            exclude_controls: Whether to exclude control samples
        
        Returns:
            DataFrame of top performers sorted by activity_score (descending)
        """
        result = df.copy()
        
        if exclude_controls:
            result = result[result['is_control'] == False]
        
        # Remove rows with NaN activity scores
        result = result.dropna(subset=['activity_score'])
        
        # Sort by activity score descending
        result = result.sort_values('activity_score', ascending=False)
        
        return result.head(n)
    
    def get_generation_statistics(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate statistics per generation for visualization
        
        Args:
            df: DataFrame with generation and activity_score columns
        
        Returns:
            DataFrame with generation statistics
        """
        # Exclude controls and NaN scores
        variants = df[
            (df['is_control'] == False) & 
            (df['activity_score'].notna())
        ].copy()
        
        if variants.empty:
            return pd.DataFrame()
        
        stats = variants.groupby('generation').agg(
            count=('activity_score', 'size'),
            mean_activity=('activity_score', 'mean'),
            median_activity=('activity_score', 'median'),
            min_activity=('activity_score', 'min'),
            max_activity=('activity_score', 'max'),
            std_activity=('activity_score', 'std'),
            q25_activity=('activity_score', lambda x: x.quantile(0.25)),
            q75_activity=('activity_score', lambda x: x.quantile(0.75))
        ).reset_index()
        
        # Fill NaN std with 0 (happens when only 1 sample in generation)
        stats['std_activity'] = stats['std_activity'].fillna(0)
        
        return stats


# Singleton instance
activity_calculator = ActivityScoreCalculator()
