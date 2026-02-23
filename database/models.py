from dataclasses import dataclass
from typing import Optional


@dataclass
class MutationInfo:
    """Stores information about a single codon-level mutation."""
    position: int
    wt_codon: str
    mut_codon: str
    wt_aa: str
    mut_aa: str
    mutation_type: str  # 'synonymous' or 'non-synonymous'

    def __str__(self):
        return f"{self.wt_aa}{self.position}{self.mut_aa} ({self.mutation_type})"


@dataclass
class StagingResult:
    """
    Mirrors the staging_results table produced by the plasmid validation pipeline.
    Describes where (and how confidently) the gene was located in the plasmid.

    Fields
    ------
    id                : primary key from staging_results
    experiment_id     : links back to the experiment
    is_valid          : True if the validation pipeline accepted this result
    match_type        : e.g. 'exact', 'high_identity', 'partial'
    identity          : alignment identity as a fraction 0-1
    coverage          : alignment coverage as a fraction 0-1
    strand            : '+' (forward) or '-' (reverse)
    frame             : reading frame offset — 0, 1, or 2
    start_nt          : 0-based start position on the canonical plasmid
    end_nt_exclusive  : exclusive end position
    wraps_origin      : True if the gene spans the plasmid origin
    notes             : optional free-text from the validation pipeline
    created_at        : ISO-8601 timestamp string
    """
    id: int
    experiment_id: int
    is_valid: bool
    match_type: str
    identity: float
    coverage: float
    strand: str
    frame: int
    start_nt: int
    end_nt_exclusive: int
    wraps_origin: bool
    notes: Optional[str]
    created_at: str
