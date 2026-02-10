import sqlite3
from typing import List, Tuple, Optional
from dataclasses import dataclass
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


@dataclass
class MutationInfo:
    """Data class to store mutation information"""
    position: int
    wt_codon: str
    mut_codon: str
    wt_aa: str
    mut_aa: str
    mutation_type: str  # 'synonymous' or 'non-synonymous'
    
    def __str__(self):
        return f"{self.wt_aa}{self.position}{self.mut_aa} ({self.mutation_type})"


class SequenceAnalyzer:
    """Sequence analysis and mutation detection using Smith-Waterman alignment"""
    
    # Standard genetic code
    GENETIC_CODE = {
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    
    def __init__(self, wt_protein_seq: str, wt_plasmid_seq: str):
        """Initialize with wild-type sequences"""
        self.wt_protein_seq = wt_protein_seq.upper().strip()
        self.wt_plasmid_seq = wt_plasmid_seq.upper().strip()
        self.wt_gene_seq = None
        self.wt_gene_start = None
        
        self._locate_wt_gene()
    
    def _smith_waterman(self, seq1: str, seq2: str, match=2, mismatch=-1, gap=-1) -> Tuple[int, int, int]:
        """
        Smith-Waterman local alignment algorithm
        Returns: (score, start_pos_in_seq1, end_pos_in_seq1)
        """
        m, n = len(seq1), len(seq2)
        
        # Initialize scoring matrix
        H = [[0] * (n + 1) for _ in range(m + 1)]
        max_score = 0
        max_pos = (0, 0)
        
        # Fill matrix
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                match_score = H[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
                delete = H[i-1][j] + gap
                insert = H[i][j-1] + gap
                H[i][j] = max(0, match_score, delete, insert)
                
                if H[i][j] > max_score:
                    max_score = H[i][j]
                    max_pos = (i, j)
        
        # Traceback to find alignment boundaries
        i, j = max_pos
        while i > 0 and j > 0 and H[i][j] > 0:
            i -= 1
            j -= 1
        
        start_pos = i
        end_pos = max_pos[0]
        
        return max_score, start_pos, end_pos
    
    def _translate(self, dna_seq: str) -> str:
        """Translate DNA to protein sequence"""
        protein = []
        for i in range(0, len(dna_seq) - 2, 3):
            codon = dna_seq[i:i+3]
            aa = self.GENETIC_CODE.get(codon, 'X')
            if aa == '*':
                break
            protein.append(aa)
        return ''.join(protein)
    
    def _reverse_complement(self, seq: str) -> str:
        """Generate reverse complement of DNA sequence"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return ''.join(complement.get(base, base) for base in reversed(seq))
    
    def _locate_wt_gene(self):
        """Locate gene using Smith-Waterman alignment on translated sequences"""
        logger.info("Locating WT gene using Smith-Waterman alignment...")
        
        # Extend plasmid for circular DNA
        extended_plasmid = self.wt_plasmid_seq + self.wt_plasmid_seq
        best_score = 0
        best_match = None
        
        # Try both strands
        for strand_seq in [extended_plasmid, self._reverse_complement(extended_plasmid)]:
            # Try all three reading frames
            for frame in range(3):
                # Translate entire sequence in this frame
                translated = self._translate(strand_seq[frame:])
                
                # Use Smith-Waterman to find best alignment with WT protein
                score, start, end = self._smith_waterman(translated, self.wt_protein_seq)
                
                if score > best_score:
                    best_score = score
                    # Calculate DNA coordinates
                    dna_start = frame + (start * 3)
                    dna_end = frame + (end * 3)
                    gene_seq = strand_seq[dna_start:dna_end]
                    
                    best_match = {
                        'score': score,
                        'gene_seq': gene_seq,
                        'start': dna_start % len(self.wt_plasmid_seq),
                        'protein': self._translate(gene_seq)
                    }
        
        if not best_match or best_match['score'] < len(self.wt_protein_seq) * 1.5:
            raise ValueError("Could not locate WT gene in plasmid!")
        
        self.wt_gene_seq = best_match['gene_seq']
        self.wt_gene_start = best_match['start']
        
        logger.info(f"WT gene found: score={best_match['score']}, start={self.wt_gene_start}, "
                   f"length={len(self.wt_gene_seq)}")
    
    def extract_gene_from_variant(self, variant_plasmid_seq: str) -> str:
        """Extract gene from variant plasmid using coordinates from WT"""
        variant_plasmid_seq = variant_plasmid_seq.upper().strip()
        gene_len = len(self.wt_gene_seq)
        
        # Extract gene (handle circular DNA)
        if self.wt_gene_start + gene_len <= len(variant_plasmid_seq):
            gene_seq = variant_plasmid_seq[self.wt_gene_start:self.wt_gene_start + gene_len]
        else:
            # Gene wraps around
            gene_seq = (variant_plasmid_seq[self.wt_gene_start:] + 
                       variant_plasmid_seq[:gene_len - (len(variant_plasmid_seq) - self.wt_gene_start)])
        
        return gene_seq
    
    def identify_mutations(self, variant_gene: str) -> List[MutationInfo]:
        """Identify mutations comparing variant to WT"""
        mutations = []
        
        if len(variant_gene) != len(self.wt_gene_seq):
            logger.warning(f"Length mismatch: variant={len(variant_gene)}, WT={len(self.wt_gene_seq)}")
            return mutations
        
        # Compare codon by codon
        for i in range(0, len(self.wt_gene_seq) - 2, 3):
            wt_codon = self.wt_gene_seq[i:i+3]
            mut_codon = variant_gene[i:i+3]
            
            if wt_codon != mut_codon:
                wt_aa = self.GENETIC_CODE.get(wt_codon, 'X')
                mut_aa = self.GENETIC_CODE.get(mut_codon, 'X')
                
                mutation = MutationInfo(
                    position=(i // 3) + 1,
                    wt_codon=wt_codon,
                    mut_codon=mut_codon,
                    wt_aa=wt_aa,
                    mut_aa=mut_aa,
                    mutation_type='synonymous' if wt_aa == mut_aa else 'non-synonymous'
                )
                mutations.append(mutation)
        
        return mutations
    
    def analyze_variant(self, variant_plasmid: str) -> dict:
        """
        Complete analysis of a variant
        Returns dict with gene_seq, protein_seq, and mutations
        """
        gene_seq = self.extract_gene_from_variant(variant_plasmid)
        protein_seq = self._translate(gene_seq)
        mutations = self.identify_mutations(gene_seq)
        
        num_synonymous = sum(1 for m in mutations if m.mutation_type == 'synonymous')
        num_nonsynonymous = len(mutations) - num_synonymous
        
        logger.info(f"Analysis complete: {len(mutations)} mutations "
                   f"({num_synonymous} syn, {num_nonsynonymous} non-syn)")
        
        return {
            'gene_sequence': gene_seq,
            'protein_sequence': protein_seq,
            'mutations': mutations,
            'num_mutations': len(mutations),
            'num_synonymous': num_synonymous,
            'num_nonsynonymous': num_nonsynonymous
        }


class SequenceDatabase:
    """Database handler for storing sequence analysis results"""
    
    def __init__(self, db_path: str = "directed_evolution.db"):
        self.db_path = db_path
        self.conn = sqlite3.connect(self.db_path)
        self.cursor = self.conn.cursor()
        self._create_tables()
        logger.info(f"Connected to database: {self.db_path}")
    
    def _create_tables(self):
        """Create database tables"""
        self.cursor.execute("""
            CREATE TABLE IF NOT EXISTS variants (
                variant_id INTEGER PRIMARY KEY,
                gene_sequence TEXT,
                protein_sequence TEXT,
                num_mutations INTEGER,
                num_synonymous INTEGER,
                num_nonsynonymous INTEGER
            )
        """)
        
        self.cursor.execute("""
            CREATE TABLE IF NOT EXISTS mutations (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                variant_id INTEGER,
                position INTEGER,
                wt_codon TEXT,
                mut_codon TEXT,
                wt_aa TEXT,
                mut_aa TEXT,
                mutation_type TEXT,
                FOREIGN KEY (variant_id) REFERENCES variants(variant_id)
            )
        """)
        
        self.conn.commit()
    
    def store_variant(self, variant_id: int, analysis: dict):
        """Store variant analysis in database"""
        self.cursor.execute("""
            INSERT OR REPLACE INTO variants 
            (variant_id, gene_sequence, protein_sequence, num_mutations, 
             num_synonymous, num_nonsynonymous)
            VALUES (?, ?, ?, ?, ?, ?)
        """, (
            variant_id,
            analysis['gene_sequence'],
            analysis['protein_sequence'],
            analysis['num_mutations'],
            analysis['num_synonymous'],
            analysis['num_nonsynonymous']
        ))
        
        # Store mutations
        for mutation in analysis['mutations']:
            self.cursor.execute("""
                INSERT INTO mutations 
                (variant_id, position, wt_codon, mut_codon, wt_aa, mut_aa, mutation_type)
                VALUES (?, ?, ?, ?, ?, ?, ?)
            """, (
                variant_id,
                mutation.position,
                mutation.wt_codon,
                mutation.mut_codon,
                mutation.wt_aa,
                mutation.mut_aa,
                mutation.mutation_type
            ))
        
        self.conn.commit()
        logger.info(f"Stored analysis for variant {variant_id}")
    
    def get_variant(self, variant_id: int) -> Optional[dict]:
        """Retrieve variant analysis from database"""
        self.cursor.execute("SELECT * FROM variants WHERE variant_id = ?", (variant_id,))
        row = self.cursor.fetchone()
        
        if not row:
            return None
        
        # Get mutations
        self.cursor.execute("SELECT * FROM mutations WHERE variant_id = ?", (variant_id,))
        mutations = self.cursor.fetchall()
        
        return {
            'variant_id': row[0],
            'gene_sequence': row[1],
            'protein_sequence': row[2],
            'num_mutations': row[3],
            'num_synonymous': row[4],
            'num_nonsynonymous': row[5],
            'mutations': mutations
        }
    
    def close(self):
        """Close database connection"""
        self.conn.close()
        logger.info("Database connection closed")


