# Part C — Staging: UniProt Integration & Plasmid Validation

## Overview
This module stages a directed evolution experiment by linking a wild‑type (WT) protein to a plasmid that encodes it. The user provides a UniProt accession and a plasmid FASTA. The system retrieves the canonical WT protein sequence and feature annotations from UniProt, then validates that the WT protein is encoded in the plasmid sequence. This ensures the experiment is grounded on a correct genotype–phenotype link.

## UniProt Integration
We use the UniProt REST API (`/uniprotkb/{accession}.fasta` and `/uniprotkb/{accession}.json`) to obtain:
- The protein sequence (FASTA) for the WT enzyme.
- The feature table (domains, motifs, functional regions) for downstream reporting.

## Plasmid Validation Strategy
Plasmids are circular; therefore, the gene may wrap around the end of the linear FASTA. To handle this, the plasmid DNA is doubled (`dna2 = dna + dna`) and translated in all six reading frames. WT detection proceeds via a tiered strategy:

1. Exact match with X‑wildcards
Searches for an exact WT protein substring in any frame (allowing ambiguous translation as `X`).

2. Fuzzy match
If no exact match, we compute WT‑length window identity in each frame and accept windows above a threshold (default 95%).

3. Local alignment (Smith–Waterman)
If still no hit, we perform a lightweight local alignment to support small indels, accepting matches that meet minimum identity and coverage thresholds.

This makes the approach generalizable to any UniProt protein and any plasmid, while still robust to sequencing ambiguity, indels, and circularity.

## Outputs
The staging service returns:
- WT protein sequence
- UniProt feature table (if requested)
- Validation result (location, frame, strand, match type, identity/coverage)
- Error message if UniProt fails or validation fails

## Web Integration
We expose staging in two ways:
- UI: `/staging/` form (accession + plasmid upload + validation results)
- API: `POST /staging/api/staging` (JSON in/out)
