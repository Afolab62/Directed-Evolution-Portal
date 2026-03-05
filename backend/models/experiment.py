from datetime import datetime
from sqlalchemy import Column, String, DateTime, Text, ForeignKey, Float, Boolean, Integer
from sqlalchemy.dialects.postgresql import UUID, JSONB
from sqlalchemy.orm import relationship
from database import Base
import uuid
import math


def safe_float(value):
    """Convert NaN/Infinity to None for JSON serialization
    Pandas calculations can produce
    NaN or Inf values that are valid in Python/numpy but will cause JSON
    serialization to fail."""

    if value is None:
        return None
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    return value


class Experiment(Base):
    """
    Represents a single directed evolution experiment.
    Created when a user provides a UniProt accession and a plasmid sequence.
    The plasmid is validated against the WT protein before the experiment is saved.
    """
    __tablename__ = 'experiments'

    # UUID primary key — avoids sequential integer IDs being guessable in API routes

    id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    user_id = Column(UUID(as_uuid=True), ForeignKey('users.id', ondelete='CASCADE'), nullable=False, index=True)     # Cascade delete ensures all variants are removed if the owning user is deleted

    name = Column(String(255), nullable=False)
    
    # UniProt data — fetched from UniProt REST API at experiment creation time

    protein_accession = Column(String(50), nullable=False)
    wt_protein_sequence = Column(Text, nullable=False)
    protein_features = Column(JSONB, nullable=True)      # Stores domain annotations from UniProt as JSON for use in the frontend

    # Plasmid data
    plasmid_name = Column(String(255), nullable=True)
    plasmid_sequence = Column(Text, nullable=False)
    
    # Validation results
    validation_status = Column(String(20), nullable=False, default='invalid')  # valid, invalid
    validation_message = Column(Text, nullable=True)
    validation_data = Column(JSONB, nullable=True)  # Full validation result from plasmid_validation
    
    # Analysis status for sequence/mutation analysis
    analysis_status = Column(String(20), nullable=False, default='not_started')  # not_started, analyzing, completed, failed
    analysis_message = Column(Text, nullable=True)
    
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)
    
    def to_dict(self, include_sequences=False):
        """Convert experiment to dictionary"""
        # Extract protein name from features if available
        protein_name = None
        if self.protein_features and isinstance(self.protein_features, dict):
            protein_name = self.protein_features.get('name')
        
        data = {
            'id': str(self.id),
            'userId': str(self.user_id),
            'name': self.name,
            'proteinAccession': self.protein_accession,
            'proteinName': protein_name,
            'plasmidName': self.plasmid_name,
            'validationStatus': self.validation_status,
            'validationMessage': self.validation_message,
            'analysisStatus': self.analysis_status,
            'analysisMessage': self.analysis_message,
            'createdAt': self.created_at.isoformat(),
            'updatedAt': self.updated_at.isoformat(),
        }
        
        if include_sequences:
            data['wtProteinSequence'] = self.wt_protein_sequence
            data['plasmidSequence'] = self.plasmid_sequence
            data['proteinFeatures'] = self.protein_features
            data['validationData'] = self.validation_data
        
        return data


class VariantData(Base):
    """
    Stores a single variant (or control) row from an uploaded TSV/JSON data file.
    Each variant belongs to one experiment and one generation of directed evolution.
    Controls (is_control=True) are Gen 0 wild-type references used to normalise
    activity scores across the rest of the variants.
    """
    __tablename__ = 'variant_data'
    
    id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    experiment_id = Column(UUID(as_uuid=True), ForeignKey('experiments.id', ondelete='CASCADE'), nullable=False, index=True)
    
    # Essential fields from TSV
    plasmid_variant_index = Column(Float, nullable=False, index=True)
    parent_plasmid_variant = Column(Float, nullable=True)
    generation = Column(Integer, nullable=False, index=True)
    assembled_dna_sequence = Column(Text, nullable=False)
    dna_yield = Column(Float, nullable=False)  # DNA_Quantification_fg
    protein_yield = Column(Float, nullable=False)  # Protein_Quantification_pg
    is_control = Column(Boolean, nullable=False, default=False, index=True)
    
    # Computed fields
    protein_sequence = Column(Text, nullable=True)  # Translated from DNA
    activity_score = Column(Float, nullable=True, index=True)
    
    # QC status
    qc_status = Column(String(20), nullable=False, default='pending')  # passed, failed
    qc_message = Column(Text, nullable=True)
    
    # Any extra columns from the TSV that don't map to a fixed schema field
    extra_metadata = Column(JSONB, nullable=True)
    
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    
    # One variant has many mutations. cascade="all, delete-orphan" ensures
    # mutations are deleted when the variant is deleted (no orphaned rows).
    # lazy="select" means mutations are only fetched from DB when accessed.
    # Relationships
    
    mutations = relationship("Mutation", back_populates="variant", cascade="all, delete-orphan", lazy="select")
    
    def to_dict(self, include_sequences=False, include_mutations=False):
        """Convert variant to dictionary"""
        data = {
            'id': str(self.id),
            'experimentId': str(self.experiment_id),
            'plasmidVariantIndex': safe_float(self.plasmid_variant_index),
            'parentPlasmidVariant': safe_float(self.parent_plasmid_variant),
            'generation': self.generation,
            'dnaYield': safe_float(self.dna_yield),
            'proteinYield': safe_float(self.protein_yield),
            'activityScore': safe_float(self.activity_score),
            'isControl': self.is_control,
            'qcStatus': self.qc_status,
            'qcMessage': self.qc_message,
            'metadata': self.extra_metadata or {},
            'createdAt': self.created_at.isoformat(),
            # Always include whether sequence analysis has been run
            'proteinSequence': self.protein_sequence,
        }

        if include_mutations:
            # Explicitly requested — load and serialise mutations
            data['mutations'] = [m.to_dict() for m in self.mutations] if self.mutations else []
            data['mutationCount'] = len(data['mutations'])
        else:
            # Avoid firing a lazy-load query just for the count.
            # SQLAlchemy stores loaded relationship collections in __dict__.
            if 'mutations' in self.__dict__:
                data['mutationCount'] = len(self.__dict__['mutations'])
            else:
                data['mutationCount'] = 0

        if include_sequences:
            data['assembledDNASequence'] = self.assembled_dna_sequence
            # proteinSequence already included by default above
        
        return data


class Mutation(Base):
    """Model for storing mutations found in variants"""
    __tablename__ = 'mutations'
    
    id = Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    variant_id = Column(UUID(as_uuid=True), ForeignKey('variant_data.id', ondelete='CASCADE'), nullable=False, index=True)
    
    position = Column(Integer, nullable=False)
    wild_type = Column(String(1), nullable=False)  # Single amino acid (WT)
    mutant = Column(String(1), nullable=False)  # Single amino acid (mutant)
    wt_codon = Column(String(3), nullable=True)  # Wild-type codon
    mut_codon = Column(String(3), nullable=True)  # Mutant codon
    mut_aa = Column(String(1), nullable=True)  # Mutant amino acid (for clarity, same as mutant)
    mutation_type = Column(String(20), nullable=False)  # synonymous, non-synonymous
    generation_introduced = Column(Integer, nullable=False)  # When mutation first appeared
    
    # Relationships
    variant = relationship("VariantData", back_populates="mutations")
    
    def to_dict(self):
        """Convert mutation to dictionary"""
        return {
            'position': self.position,
            'wildType': self.wild_type,
            'mutant': self.mutant,
            'wtCodon': self.wt_codon,
            'mutCodon': self.mut_codon,
            'mutAa': self.mut_aa,
            'type': self.mutation_type,
            'generation': self.generation_introduced,
        }
