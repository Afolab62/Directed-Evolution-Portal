# Part C — Plasmid Validation (Staging)

This module stages a directed evolution experiment by validating that an uploaded plasmid FASTA encodes the selected UniProt WT protein. It is generalizable to any UniProt accession and plasmid FASTA and handles circular DNA.

---

## Requirements
- Python 3.11+  
- Flask + dependencies installed via `requirements.txt`

---

## What This Does
- Fetches WT protein FASTA and feature table from UniProt
- Validates that the plasmid encodes the WT protein
- Handles circular plasmids (wraparound)
- Provides UI + JSON API for the same workflow
- Returns clear diagnostics on failure

---

## How Validation Works 
1. Circular handling: plasmid is doubled to detect wraparound genes  
2. Six‑frame translation: check all strands/frames  
3. Tiered matching: exact → fuzzy → Smith–Waterman (indel‑capable)  

---

## Demo — Success Case
Accession: O34996  
FASTA: data/pET-28a_BSU_DNA_Pol_I_WT.fa  

Steps:
1) Run the server  
2) Open /staging/  
3) Click Load Example FASTA  
4) Click Validate  

Expected: is_valid: true, match_type: exact

---

## Demo — Fuzzy Success Case (Variant)
Accession: P61875  
FASTA: use the matching plasmid FASTA you tested (fuzzy ~0.994 identity)  

Expected: is_valid: true, match_type: fuzzy

---

## Demo — Error Message Case
Accession: BADACC  
FASTA (paste):  
    >test  
    ACGTACGTACGTACGTACGTACGTACGTACGT  

Expected: UniProt error returned cleanly.

---

## Demo — Stress Test (Large FASTA)
Accession: P61875  
FASTA: pXO1.fasta (~181 kb)  

Expected: clean failure diagnostics + alignment skipped due to size guard.

---

## Run the App
    cd /Users/munaberhe/Directed-Evolution-Portal  
    source .venv/bin/activate  
    flask --app app:create_app run  

Open in browser:  
    http://127.0.0.1:5000/staging/

---

## API Usage
    curl -X POST http://127.0.0.1:5000/staging/api/staging \
      -H "Content-Type: application/json" \
      -d '{"accession":"O34996","plasmid_fasta":">p\nACGT...","fetch_features":false}'

---

## Output Fields (Key)
- is_valid — true/false  
- match_type — exact / fuzzy / align / none  
- identity, coverage  
- frame, strand  
- start_nt, end_nt_exclusive  
- wraps_origin  

---

## Troubleshooting
- Internal Server Error (Flask‑Login): add a temporary `user_loader` in `app/__init__.py` for demo mode.  
- UniProt errors: invalid accession or network issue — try a known accession (e.g., O34996).  
- Large FASTA takes long: alignment step is skipped by size guard for safety.

---

## Limitations
- Smith–Waterman alignment is O(nm); size guard skips alignment for very large inputs.  
- Validation requires the plasmid to actually encode the chosen UniProt protein.  
- UniProt feature table is stored as‑is without further curation.

---

## References
- UniProt Consortium. *UniProt: the Universal Protein Knowledgebase in 2025*. **Nucleic Acids Research** (2025). DOI: 10.1093/nar/gkaf005  
- Smillie C. et al. *Mobility of plasmids*. **Microbiology and Molecular Biology Reviews** (2010). DOI: 10.1128/MMBR.00020-10  
- Smith TF, Waterman MS. *Identification of common molecular subsequences*. **Journal of Molecular Biology** (1981). DOI: 10.1016/0022-2836(81)90087-5  
- Rice P, Longden I, Bleasby A. *EMBOSS: The European Molecular Biology Open Software Suite*. **Trends in Genetics** (2000). DOI: 10.1016/S0168-9525(00)02024-2  
