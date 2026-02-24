"""
database/models.py
------------------
Shared dataclasses used across the entire project.
These are plain Python objects — no DB logic here.
"""

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class StagingResult:
    """
    Represents a validated staging result row from the staging_results table.
    Produced by the plasmid validation pipeline and consumed by SequenceAnalyzer.
    """
    id: int
    experiment_id: int
    is_valid: bool
    match_type: str
    identity: float
    coverage: float
    strand: str           # '+' or '-'
    frame: int            # 0, 1, or 2
    start_nt: int
    end_nt_exclusive: int
    wraps_origin: bool
    notes: Optional[str] = None


@dataclass
class MutationInfo:
    """
    Represents a single codon-level mutation identified by SequenceAnalyzer.
    """
    position: int          # 1-based amino acid position
    wt_codon: str
    mut_codon: str
    wt_aa: str
    mut_aa: str
    mutation_type: str     # 'synonymous' or 'non-synonymous'
