import sqlite3
import logging
from typing import List, Optional

try:
    from .models import MutationInfo, StagingResult
except ImportError:  # pragma: no cover
    from models import MutationInfo, StagingResult


logger = logging.getLogger(__name__)


class SequenceDatabase:
    """
    Database handler for storing sequence analysis results.

    Also provides helpers for reading from the staging_results table
    managed by the plasmid validation pipeline.
    """

    def __init__(self, db_path: str = "database.db"):
        self.db_path = db_path
        self.conn = sqlite3.connect(self.db_path)
        self.conn.row_factory = sqlite3.Row  # column access by name
        self.cursor = self.conn.cursor()
        self._create_tables()
        logger.info(f"Connected to database: {self.db_path}")

    # ------------------------------------------------------------------
    # Schema
    # ------------------------------------------------------------------

    def _create_tables(self) -> None:
        """Create analysis tables. (staging_results is owned by the validation pipeline.)"""
        self.cursor.executescript("""
            CREATE TABLE IF NOT EXISTS variants (
                variant_id          INTEGER PRIMARY KEY,
                experiment_id       INTEGER,
                staging_result_id   INTEGER REFERENCES staging_results(id),
                gene_sequence       TEXT,
                protein_sequence    TEXT,
                num_mutations       INTEGER,
                num_synonymous      INTEGER,
                num_nonsynonymous   INTEGER
            );

            CREATE TABLE IF NOT EXISTS mutations (
                id              INTEGER PRIMARY KEY AUTOINCREMENT,
                variant_id      INTEGER REFERENCES variants(variant_id),
                position        INTEGER,
                wt_codon        TEXT,
                mut_codon       TEXT,
                wt_aa           TEXT,
                mut_aa          TEXT,
                mutation_type   TEXT
            );
        """)
        self.conn.commit()

    # ------------------------------------------------------------------
    # staging_results — read side
    # ------------------------------------------------------------------

    def get_staging_result(
        self,
        experiment_id: int,
        valid_only: bool = True,
    ) -> Optional[StagingResult]:
        """
        Fetch the most recent staging result for a given experiment.

        Parameters
        ----------
        experiment_id : int
            The experiment whose plasmid was validated.
        valid_only : bool
            If True (default), only return rows where is_valid = 1.
            Set to False to retrieve failed/flagged results for inspection.
        """
        query = """
            SELECT *
            FROM   staging_results
            WHERE  experiment_id = ?
            {}
            ORDER  BY created_at DESC
            LIMIT  1
        """.format("AND is_valid = 1" if valid_only else "")

        self.cursor.execute(query, (experiment_id,))
        row = self.cursor.fetchone()

        if row is None:
            return None

        return self._row_to_staging_result(row)

    def get_all_staging_results(self, experiment_id: int) -> List[StagingResult]:
        """Return all staging results for an experiment (useful for debugging)."""
        self.cursor.execute(
            "SELECT * FROM staging_results WHERE experiment_id = ? ORDER BY created_at DESC",
            (experiment_id,)
        )
        return [self._row_to_staging_result(r) for r in self.cursor.fetchall()]

    @staticmethod
    def _row_to_staging_result(row) -> StagingResult:
        return StagingResult(
            id=row['id'],
            experiment_id=row['experiment_id'],
            is_valid=bool(row['is_valid']),
            match_type=row['match_type'],
            identity=row['identity'],
            coverage=row['coverage'],
            strand=row['strand'],
            frame=row['frame'],
            start_nt=row['start_nt'],
            end_nt_exclusive=row['end_nt_exclusive'],
            wraps_origin=bool(row['wraps_origin']),
            notes=row['notes'],
            created_at=row['created_at'],
        )

    # ------------------------------------------------------------------
    # variants / mutations — write / read
    # ------------------------------------------------------------------

    def store_variant(
        self,
        variant_id: int,
        analysis: dict,
        experiment_id: Optional[int] = None,
        staging_result_id: Optional[int] = None,
    ) -> None:
        """Store variant analysis, optionally linking back to its staging result."""
        self.cursor.execute("""
            INSERT OR REPLACE INTO variants
                (variant_id, experiment_id, staging_result_id,
                 gene_sequence, protein_sequence,
                 num_mutations, num_synonymous, num_nonsynonymous)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            variant_id,
            experiment_id,
            staging_result_id,
            analysis['gene_sequence'],
            analysis['protein_sequence'],
            analysis['num_mutations'],
            analysis['num_synonymous'],
            analysis['num_nonsynonymous'],
        ))

        for mutation in analysis['mutations']:
            self.cursor.execute("""
                INSERT INTO mutations
                    (variant_id, position, wt_codon, mut_codon,
                     wt_aa, mut_aa, mutation_type)
                VALUES (?, ?, ?, ?, ?, ?, ?)
            """, (
                variant_id,
                mutation.position,
                mutation.wt_codon,
                mutation.mut_codon,
                mutation.wt_aa,
                mutation.mut_aa,
                mutation.mutation_type,
            ))

        self.conn.commit()
        logger.info(f"Stored analysis for variant {variant_id}")

    def get_variant(self, variant_id: int) -> Optional[dict]:
        """Retrieve variant analysis from the database."""
        self.cursor.execute("SELECT * FROM variants WHERE variant_id = ?", (variant_id,))
        row = self.cursor.fetchone()
        if not row:
            return None

        self.cursor.execute("SELECT * FROM mutations WHERE variant_id = ?", (variant_id,))
        mutations = self.cursor.fetchall()

        return {
            'variant_id':        row['variant_id'],
            'experiment_id':     row['experiment_id'],
            'staging_result_id': row['staging_result_id'],
            'gene_sequence':     row['gene_sequence'],
            'protein_sequence':  row['protein_sequence'],
            'num_mutations':     row['num_mutations'],
            'num_synonymous':    row['num_synonymous'],
            'num_nonsynonymous': row['num_nonsynonymous'],
            'mutations':         [dict(m) for m in mutations],
        }

    def close(self) -> None:
        """Close database connection."""
        self.conn.close()
        logger.info("Database connection closed")
