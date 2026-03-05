// Core types for the Direct Evolution Monitoring Portal

export interface User {
  id: string;
  email: string;
  passwordHash: string;
  createdAt: Date;
}

export interface UniProtProtein {
  accession: string;
  name: string;
  organism: string;
  sequence: string;
  length: number;
  features: UniProtFeature[];
  // Enriched fields returned by the detailed endpoint
  gene_names?: string[];
  keywords?: string[];
  taxonomy_lineage?: string[];
  pdb_ids?: string[];
  alphafold_id?: string | null;
  interpro_ids?: string[];
  pfam_ids?: string[];
  kegg_ids?: string[];
  go_terms?: { id: string; term: string; aspect: string }[];
  crossrefs?: Record<string, string[]>;
}

export interface UniProtFeature {
  type: string;
  description: string;
  location: {
    start: number;
    end: number;
  };
}

export interface Experiment {
  id: string;
  userId: string;
  name: string;
  proteinAccession: string;
  proteinName?: string | null;
  protein: UniProtProtein | null;
  wtProteinSequence?: string;
  proteinFeatures?: any;
  plasmidSequence: string;
  plasmidName: string;
  validationStatus: "pending" | "valid" | "invalid";
  validationMessage: string;
  analysisStatus: "not_started" | "analyzing" | "completed" | "failed";
  analysisMessage?: string | null;
  createdAt: Date;
  updatedAt: Date;
}

export interface VariantData {
  id: string;
  experimentId: string;
  plasmidVariantIndex: number;
  parentPlasmidVariant: number | null;
  generation: number;
  assembledDNASequence?: string;
  proteinSequence: string | null;
  dnaYield: number;
  proteinYield: number;
  activityScore: number;
  isControl: boolean;
  mutations?: Mutation[]; // Optional: only available after sequence analysis
  mutationCount?: number; // Always available after analysis (even without full mutation details)
  qcStatus: "passed" | "failed";
  qcMessage: string;
  metadata: Record<string, unknown>;
}

export interface Mutation {
  position: number;
  wildType: string;
  mutant: string;
  wtCodon?: string; // DNA codon for wild-type
  mutCodon?: string; // DNA codon for mutant
  mutAa?: string; // Mutant amino acid (alias for backward compatibility)
  type: "synonymous" | "non-synonymous";
  generation: number;
}

export interface ParseResult {
  success: boolean;
  variants: VariantData[];
  errors: ParseError[];
  warnings: string[];
}

export interface ParseError {
  row: number;
  field: string;
  message: string;
}

export interface AnalysisResult {
  topPerformers: VariantData[];
  generationStats: GenerationStats[];
  totalVariants: number;
  passedQC: number;
  failedQC: number;
}

export interface GenerationStats {
  generation: number;
  count: number;
  meanActivity: number;
  medianActivity: number;
  minActivity: number;
  maxActivity: number;
  stdDev: number;
}
