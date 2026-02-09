# Sequence Analysis Module - Section 4.a
## DNA to Protein Processing for Directed Evolution


---

## Files Included

1. **`sequence_analyzer.py`** - Core analysis module
   - `SequenceAnalyzer` class for DNA→Protein translation and mutation detection
   - `SequenceDatabase` class for storing results in SQLite database
   - Handles circular plasmid DNA sequences

#STEPS:

### Step 1: Prepare Data



### Step 2: Run the Analysis

```bash
python process_data.py
```

This will:
1. Load the WT sequences
2. Parse experimental data
3. Analyze each variant (extract gene, translate, find mutations)
4. Store everything in SQLite database
5. Generate processing summary

### Step 3: Generate Reports

```bash
python database_queries.py
```

This generates:
- `all_variants.csv` - Complete variant information
- `all_mutations.csv` - All mutations across all variants
- `generation_summary.csv` - Statistics per generation
- `mutation_frequencies.csv` - Mutation hotspot analysis

---

##  Database Schema

### Table: `variant_sequences`

| Column | Type | Description |
|--------|------|-------------|
| variant_index | INTEGER | Unique identifier (Primary Key) |
| parent_variant | INTEGER | Parent variant index |
| generation | INTEGER | Directed evolution generation |
| dna_sequence | TEXT | Full plasmid DNA sequence |
| gene_sequence | TEXT | Extracted gene sequence |
| protein_sequence | TEXT | Translated protein sequence |
| num_mutations | INTEGER | Total mutations vs WT |
| num_synonymous | INTEGER | Synonymous mutations |
| num_nonsynonymous | INTEGER | Non-synonymous mutations |

### Table: `mutations`

| Column | Type | Description |
|--------|------|-------------|
| mutation_id | INTEGER | Auto-increment ID (Primary Key) |
| variant_index | INTEGER | Foreign key to variant_sequences |
| position | INTEGER | Amino acid position (1-indexed) |
| wt_codon | TEXT | Wild-type DNA codon |
| mut_codon | TEXT | Mutant DNA codon |
| wt_aa | TEXT | Wild-type amino acid |
| mut_aa | TEXT | Mutant amino acid |
| mutation_type | TEXT | 'synonymous' or 'non-synonymous' |
| generation_introduced | INTEGER | Generation when mutation occurred |

---

## key features

### 1. Circular DNA Handling
The analyzer correctly handles circular plasmid sequences where genes may wrap around the origin:

```python
analyzer = SequenceAnalyzer(wt_protein_seq, wt_plasmid_seq)
# Automatically finds gene location accounting for circularity
```

### 2. Complete Mutation Classification
Each mutation is analyzed at both DNA and protein levels:

```python
MutationInfo:
  - position: Amino acid position
  - wt_codon/mut_codon: DNA level changes
  - wt_aa/mut_aa: Protein level changes
  - mutation_type: 'synonymous' or 'non-synonymous'
```

### 3. Lineage Tracing
Track how mutations accumulate through generations:

```python
tool = AnalysisQueryTool(db_path)
lineage = tool.get_lineage_trace(variant_index=100)
# Returns full ancestry from WT to variant 100





