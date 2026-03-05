# Directed Evolution Portal — User Help Guide

> **Quick links:** [Getting Started](#2-getting-started) · [Upload Format](#4-uploading-experimental-data) · [Activity Score](#6-activity-score-calculation) · [Sequence Analysis](#7-sequence-analysis) · [Visualisations](#8-visualisations) · [Troubleshooting](#10-troubleshooting)

---

## 1. What Is the Portal?

The Directed Evolution Portal is a web application for monitoring and analysing data produced by automated directed-evolution pipelines. It lets you:

- Register one or more **experiments**, each linked to a UniProt protein and a plasmid sequence.
- **Upload** per-variant experimental measurements (DNA yield, protein yield, generation) in TSV or JSON format.
- Automatically calculate a normalised **activity score** for every variant.
- Run **sequence analysis** to detect protein-level mutations for each variant DNA sequence.
- Explore results through interactive **visualisations**: violin distributions, mutation fingerprints, and 3-D activity landscapes.

---

## 2. Getting Started

### 2.1 Creating an Account

1. Navigate to `http://localhost:3000` (or your deployed URL).
2. Click **Register** and enter your email address and a password.  
   Passwords are hashed with bcrypt and never stored in plain text.
3. After registration you are automatically logged in and redirected to the dashboard.

### 2.2 Dashboard Overview

The dashboard shows a summary of all your experiments. Each card displays the experiment name, protein accession, validation status, and a count of uploaded variants.

---

## 3. Creating an Experiment

Click **New Experiment** from the dashboard or navigation sidebar.

**Step 1 — Protein**  
Enter a UniProt accession (e.g. `O34996`). The portal will look up the protein name, organism, and wild-type (WT) sequence from UniProt automatically. A preview of the first 100 amino acids is shown so you can confirm the correct protein was found.

**Step 2 — Plasmid**  
Upload a FASTA file (`.fasta`, `.fa`, `.fna`, or `.txt`) containing the circular plasmid DNA sequence that encodes your protein of interest. The portal validates that the uploaded plasmid actually encodes the WT protein by:

1. Translating all six reading frames (3 forward, 3 reverse-complement) of the doubled sequence (to handle origin-wrapping ORFs).
2. Searching for the WT protein as an exact amino-acid substring using fast string matching, backed by a Smith–Waterman local alignment fallback for ambiguous or low-identity sequences.
3. Checking for an ATG start codon and canonical stop codon at the called ORF boundaries.

The portal reports the validation result (_valid_ / _invalid_) along with the identified reading-frame strand, start nucleotide, and any plausibility warnings.

**Step 3 — Review and Submit**  
Confirm the experiment name, plasmid name, and validation result, then click **Create Experiment**.

---

## 4. Uploading Experimental Data

Navigate to an experiment's detail page, then select the **Upload Data** tab.

### 4.1 Supported File Formats

| Format                | Extension      |
| --------------------- | -------------- |
| Tab-separated values  | `.tsv`, `.txt` |
| JSON array of objects | `.json`        |

### 4.2 Required Columns

The following columns must be present (column headers are case-insensitive and common synonyms are accepted automatically — see §4.3):

| Canonical name           | Type    | Description                                                      |
| ------------------------ | ------- | ---------------------------------------------------------------- |
| `Plasmid_Variant_Index`  | number  | Unique identifier for the variant (matches your lab's numbering) |
| `Generation`             | integer | Directed-evolution generation (0 = wild-type controls)           |
| `Assembled_DNA_Sequence` | string  | Full-length nucleotide sequence of the assembled plasmid variant |
| `DNA_Yield`              | number  | Measured DNA yield (any consistent unit, e.g. fg)                |
| `Protein_Yield`          | number  | Measured protein yield (any consistent unit, e.g. pg)            |
| `Is_Control`             | boolean | `TRUE`/`FALSE` or `1`/`0` — marks generation-0 WT controls       |

### 4.3 Optional Columns and Accepted Synonyms

`Parent_Plasmid_Variant` (the index of the parent variant in the lineage) is optional and used to determine which generation each mutation was first introduced.

The portal automatically recognises common column-name variations, for example:

| Canonical               | Also accepted                                                 |
| ----------------------- | ------------------------------------------------------------- |
| `Plasmid_Variant_Index` | `variant_index`, `plasmid_id`, `variant_id`, `index`          |
| `Generation`            | `directed_evolution_generation`, `gen`                        |
| `DNA_Yield`             | `dna_quantification_fg`, `dna_qty_fg`, `dna_concentration_fg` |
| `Is_Control`            | `control`, `control_sample`                                   |

**Extra columns** (any columns not in the list above) are accepted without error and stored as JSON metadata on each variant row. They will not affect QC or analysis.

### 4.4 Quality Control

Every row is checked before storage:

- Required numeric fields (`Plasmid_Variant_Index`, `Generation`, `DNA_Yield`, `Protein_Yield`) must not be empty.
- `Generation` must be ≥ 0.
- `DNA_Yield` and `Protein_Yield` must be ≥ 0.
- `Is_Control` must be a recognised boolean value.
- `Assembled_DNA_Sequence` must contain only the characters A, T, C, G.

Rows that fail QC are reported in the upload response (`failedQC` count and a table of errors) but are **not discarded silently** — you will always see which rows were rejected and why. Rows passing QC are stored immediately.

> **Note:** At least one `Is_Control = TRUE` row per generation is required for activity score calculation. If a generation lacks controls, the overall median control values are used as a fallback with a warning.

---

## 5. Running Sequence Analysis

After uploading data, click **Analyse Sequences** on the experiment detail page (Upload Data tab).

This starts a **background analysis job**. The page will show a live status banner indicating the job is in progress, with an elapsed-seconds counter. The status updates every 5 seconds by polling the server — you do not need to keep the page open.

When analysis completes, a green **Analysis Complete** banner appears at the top of the page (visible on any tab) and variant mutation data becomes available in the visualisations.

### What the analysis does

1. Uses generation-0 control sequences as the WT reference plasmid (falls back to the experiment's registered plasmid if no controls are present).
2. Locates the gene of interest in the WT reference by translating all six reading frames of the doubled circular sequence, matching a protein prefix of 50 → 30 → 15 amino acids progressively.
3. **Locks the gene coordinates** (start position, length, strand) from the WT reference for all subsequent variant extractions.
4. For each variant:
   - Estimates the **circular rotation offset** — high-throughput sequencing assemblers may produce assemblies that start at a different position on the circular plasmid. The portal detects and corrects for this shift by anchored voting over short sequence windows.
   - Extracts the gene region from the variant plasmid using the corrected coordinates.
   - **Translates** the variant gene to a protein sequence using the standard genetic code.
   - Detects **mutations** by codon-by-codon comparison with the WT gene.
   - Runs **Needleman-Wunsch global alignment** of the WT and variant proteins to compute alignment-aware residue positions (important for variants with insertions or deletions near mutation sites).
5. Computes `generation_introduced` for each mutation by walking the parent–child lineage: a mutation inherited from a parent carries the parent's generation, while a newly appearing mutation is assigned the current variant's generation.

---

## 6. Activity Score Calculation

The activity score quantifies the enzymatic performance of a variant relative to the generation-matched WT controls. It is based on the principle that a more active enzyme produces **more DNA per unit of protein expressed**.

### 6.1 Formula

$$
\text{Activity Score} = \frac{\text{DNA Yield} / \text{DNA Baseline}}{\text{Protein Yield} / \text{Protein Baseline}}
$$

Where:

- **DNA Baseline** — the _median_ DNA yield of all `Is_Control = TRUE` rows in the same generation.
- **Protein Baseline** — the _median_ protein yield of all `Is_Control = TRUE` rows in the same generation.

The numerator is the **fold-change in DNA yield** relative to controls; the denominator is the **fold-change in protein expression** relative to controls. Dividing one by the other corrects for differences in how much protein the variant produces, leaving a score that reflects catalytic activity specifically.

### 6.2 Why This Formula?

| Score value | Interpretation                                                        |
| ----------- | --------------------------------------------------------------------- |
| > 1         | Variant produces more DNA per unit protein than control → more active |
| ≈ 1         | Similar activity to the WT control                                    |
| < 1         | Less active than WT control                                           |

Using fold-changes rather than raw differences makes the metric **dimensionless** and **comparable across generations** (even if absolute yield levels drift between rounds).

### 6.3 Edge Cases

- An `epsilon = 0.01` floor is applied to all yields and baselines before division, preventing division-by-zero errors from very low measurements.
- Negative yields are treated as 0 before scoring.
- Controls are assigned `activity_score = NaN` (not meaningful for controls).
- If a generation has no controls, the **overall median** of all controls is used as the baseline for that generation, with a warning logged.

---

## 7. Sequence Analysis

### 7.1 DNA → Protein Translation

The portal uses the **standard genetic code** (NCBI table 1). Translation starts from the first codon of the identified ORF and stops at the first in-frame stop codon. No alternative start codons or non-standard tables are used.

### 7.2 ORF Identification

The WT plasmid is treated as **circular** throughout. This is implemented by doubling the sequence (`dna2 = plasmid + plasmid`) before any search. Gene coordinates are then mapped back to the original plasmid (0-based, `start_nt` to `end_nt_exclusive`).

Both **sense (+) and antisense (–) strands** are searched. If the gene is on the minus strand, the extracted region is reverse-complemented before translation and mutation comparison.

### 7.3 Mutation Detection

Mutations are detected by codon-by-codon comparison of the WT and variant gene sequences:

- **Synonymous** — codon differs but encodes the same amino acid.
- **Non-synonymous** — codon encodes a different amino acid.
- Stop codons within the coding sequence are skipped (not reported as mutations).

Reporteed mutations include: 1-based protein position, WT amino acid, variant amino acid, WT codon, variant codon, mutation type, and HGVS-style change string (e.g. `A45V`).

### 7.4 Needleman-Wunsch Alignment

After codon comparison the WT and variant protein sequences are globally aligned using Needleman-Wunsch (match = +2, mismatch = −1, gap = −2). This produces an `aligned_position` for every mutation, which is the residue index on the alignment rather than the raw codon index. This is used by the 3-D mutation fingerprint to correctly overlay mutations onto AlphaFold structure residues when insertions or deletions shift the numbering.

---

## 8. Visualisations

All visualisations are on the **Analysis** page, accessible from the top navigation. Select an experiment from the dropdown to load its data.

### 8.1 Activity Score Distribution (Violin Plot)

**Tab: Overview**

A server-rendered violin plot (matplotlib) showing the full distribution of activity scores for each generation. Each violin represents the kernel-density-estimated (KDE) distribution of scores. Overlaid elements:

- Vertical line — range (min to max)
- Tick marks — minimum and maximum
- Thick horizontal bar — mean (blue)
- Thick horizontal bar — median (violet)

The summary table below the plot shows mean, median, min, max, and standard deviation per generation.

### 8.2 Mutation Fingerprint

**Tab: Mutations**

Select a variant from the dropdown. Two sub-tabs are available:

**Linear fingerprint** — A horizontal bar plot showing every mutation in the selected variant against the WT protein sequence. Mutations are represented as coloured triangles at their protein position. Colour encodes the generation in which the mutation was first introduced (lighter = earlier generation). Hover over a triangle to see the full amino-acid change and generation.

**3-D fingerprint** — The same mutations are mapped onto the AlphaFold predicted protein structure (PDB coordinates retrieved automatically). Residues are coloured by mutation count: blue = wild-type, red gradient = highly mutated. The plot is fully interactive (rotate, zoom, pan).

### 8.3 Activity Landscape

**Tab: Activity Landscape**

A 3-D scatter plot of all variants embedded in two dimensions using:

- **PCA** — Principal Component Analysis (fast, deterministic; default)
- **t-SNE** — t-distributed Stochastic Neighbour Embedding (non-linear, slower)
- **UMAP** — Uniform Manifold Approximation and Projection (requires optional `umap-learn` dependency)

The embedding is computed from binary mutation vectors (each dimension = presence/absence of a specific mutation). The z-axis and point colour both represent activity score. Per-generation animation frames allow you to step through the evolution trajectory.

---

## 9. Data Export

Currently, the portal does not provide a one-click bulk export. The following workarounds are available:

- **Raw variant data**: The Variants tab on the experiment detail page shows the first 50 variants sorted by activity score. For larger exports, query the Neon database directly.
- **Plots**: Right-click any Plotly chart → _Download plot as PNG_. The violin plot image can be right-clicked and saved directly.
- **FASTA**: On the New Experiment page, after looking up a protein, a **Download FASTA** link downloads the WT protein sequence as a `.fasta` file.

---

## 10. Troubleshooting

### "Analysis in progress" does not resolve

- Check that the Flask backend is running (`python run.py` from the `backend/` directory).
- Background analysis runs in a daemon thread. If the server was restarted mid-analysis, the `analysis_status` column may be stuck at `analyzing`. Click **Re-analyse** to restart.
- For large datasets (> 500 variants) analysis can take 2–5 minutes. The elapsed-seconds counter on the status banner confirms the job is still running.

### Upload rejected — "Could not find required columns"

- Check that your file is genuinely tab-separated (TSV). Open it in a text editor and confirm columns are separated by `\t`, not commas or spaces.
- Verify all required column names are present. Accepted synonyms are listed in §4.3.
- JSON uploads must be a valid JSON array of objects (not wrapped in a root key).

### Plasmid validation fails — "valid" but warnings about start codon

This is a non-fatal warning. The portal found the WT protein in the plasmid but the three nucleotides immediately before the match are not `ATG`. This can happen when:

- The FASTA sequence is not the full expression plasmid (e.g. only the gene insert was uploaded).
- The gene uses a non-ATG start codon (rare in expression vectors).

The portal will still create the experiment and perform sequence analysis correctly.

### Violin / fingerprint plot shows "Failed to load"

- The Flask backend must be running and accessible at `NEXT_PUBLIC_BACKEND_URL`.
- Check the browser network tab for the specific error returned by the server.
- Restart the backend and refresh the page.

### Activity scores all appear very similar (near 1.0)

This is expected when variant yields are close to control yields. It indicates the current generation has not yet diverged significantly from the WT. Scores should spread out in later generations as beneficial mutations accumulate.

---

## 11. File Format Reference

### Minimal TSV Example

```tsv
Plasmid_Variant_Index	Generation	Assembled_DNA_Sequence	DNA_Yield	Protein_Yield	Is_Control	Parent_Plasmid_Variant
0	0	ATGAAAGCAAT...	120.5	85.2	TRUE
1	1	ATGAAAGCAAT...	198.3	76.1	FALSE	0
2	1	ATGAAAGCGAT...	145.2	80.4	FALSE	0
```

### Minimal JSON Example

```json
[
  {
    "Plasmid_Variant_Index": 0,
    "Generation": 0,
    "Assembled_DNA_Sequence": "ATGAAAGCAAT...",
    "DNA_Yield": 120.5,
    "Protein_Yield": 85.2,
    "Is_Control": true
  },
  {
    "Plasmid_Variant_Index": 1,
    "Generation": 1,
    "Assembled_DNA_Sequence": "ATGAAAGCAAT...",
    "DNA_Yield": 198.3,
    "Protein_Yield": 76.1,
    "Is_Control": false,
    "Parent_Plasmid_Variant": 0
  }
]
```

---

## 12. Glossary

| Term                        | Definition                                                                                                                                                  |
| --------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Activity score**          | Normalised fold-change metric measuring catalytic efficiency relative to generation-matched controls (see §6)                                               |
| **Control**                 | A variant with `Is_Control = TRUE`; represents the unmodified WT construct in a given generation and is used as baseline for activity score normalisation   |
| **Generation**              | One round of directed-evolution mutagenesis and selection; 0 = starting WT                                                                                  |
| **ORF**                     | Open reading frame — the contiguous stretch of codons between a start and stop codon that encodes the protein                                               |
| **Non-synonymous mutation** | A codon change that alters the amino-acid sequence of the protein                                                                                           |
| **Synonymous mutation**     | A codon change that does not alter the amino-acid sequence (silent mutation)                                                                                |
| **Rotation offset**         | The circular shift in nucleotide position between how the WT plasmid and a variant plasmid were assembled; corrected automatically during sequence analysis |
| **Needleman-Wunsch**        | A global sequence alignment algorithm used to assign alignment-aware positions to mutations when indels are present                                         |
| **PCA / t-SNE / UMAP**      | Dimensionality reduction methods used to embed high-dimensional mutation vectors into 2-D for the activity landscape plot                                   |
| **UniProt accession**       | A unique identifier for a protein in the UniProt database (e.g. `O34996`)                                                                                   |
