"""
Migration script: adds Protein_Sequence column to activity_measurements.

Because the assembled plasmid sequences can start at any rotational position
(circular assembly ambiguity), fixed WT coordinates don't work for every
variant. Instead, for each variant we:
  1. Concatenate the sequence with itself (handles circularity).
  2. Find the longest ORF across all 6 reading frames.
  3. Translate and store it as the protein sequence.

This mirrors what plasmid_validation.py does when searching circular plasmids.

Usage:
    python analysis/add_protein_sequences.py
"""

import os
import sqlite3

# ── codon table (matches sequence_tools.py) ───────────────────────────────────

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate_dna(seq: str) -> str:
    aa = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        aa.append(CODON_TABLE.get(codon, "X"))
    return "".join(aa)


def reverse_complement(seq: str) -> str:
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(comp.get(b, "N") for b in reversed(seq))


def longest_orf_in_circular(seq: str, min_aa: int = 200) -> str:
    """
    Find the longest ORF (M...stop) in a circular DNA sequence by searching
    all 6 reading frames of (seq + seq). Returns the protein sequence
    (stop codon stripped), or empty string if none found above min_aa.
    """
    seq = seq.upper()
    circular = seq + seq           # doubles the sequence to capture wraparound ORFs
    rc_circular = reverse_complement(seq) + reverse_complement(seq)

    best_protein = ""
    best_len = min_aa

    for search_seq in (circular, rc_circular):
        for frame in range(3):
            protein = translate_dna(search_seq[frame:])
            i = 0
            while i < len(protein):
                if protein[i] == "M":
                    j = protein.find("*", i)
                    end = j if j != -1 else len(protein)
                    orf_len = end - i
                    if orf_len > best_len:
                        best_len = orf_len
                        best_protein = protein[i:end]
                    i = (j + 1) if j != -1 else len(protein)
                else:
                    i += 1

    return best_protein


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.join(script_dir, "..")
    db_path = os.path.join(repo_root, "instance", "app.db")

    if not os.path.exists(db_path):
        raise FileNotFoundError(f"Database not found: {db_path}")

    print("Connecting to database...")
    with sqlite3.connect(db_path) as conn:
        cur = conn.cursor()

        # add column if missing
        cur.execute('PRAGMA table_info("activity_measurements")')
        cols = [row[1] for row in cur.fetchall()]
        if "Protein_Sequence" not in cols:
            print("  Adding Protein_Sequence column...")
            cur.execute('ALTER TABLE activity_measurements ADD COLUMN Protein_Sequence TEXT')
        else:
            print("  Protein_Sequence column exists — re-populating...")

        cur.execute('SELECT rowid, Plasmid_Variant_Index, Assembled_DNA_Sequence FROM activity_measurements')
        rows = cur.fetchall()
        print(f"  Processing {len(rows)} variants...\n")

        updated, failed, short = 0, 0, 0
        for rowid, variant_id, dna_seq in rows:
            try:
                protein = longest_orf_in_circular(str(dna_seq))
                if not protein:
                    print(f"  Warning: no ORF >= 200 aa found for variant {variant_id}")
                    short += 1
                    continue
                cur.execute(
                    'UPDATE activity_measurements SET Protein_Sequence=? WHERE rowid=?',
                    (protein, rowid)
                )
                updated += 1
            except Exception as e:
                print(f"  Error: variant {variant_id} — {e}")
                failed += 1

        conn.commit()

    print(f"\nDone.")
    print(f"  Updated : {updated}")
    print(f"  No ORF  : {short}")
    print(f"  Errors  : {failed}")

    # quick sanity check
    print("\nSample protein sequences:")
    with sqlite3.connect(db_path) as conn:
        rows = conn.execute(
            'SELECT Plasmid_Variant_Index, Control, Protein_Sequence '
            'FROM activity_measurements LIMIT 5'
        ).fetchall()
        for r in rows:
            p = str(r[2]) if r[2] else "NULL"
            print(f"  Variant {r[0]} | Control={r[1]} | Len={len(p)} | {p[:50]}...")


if __name__ == "__main__":
    main()
