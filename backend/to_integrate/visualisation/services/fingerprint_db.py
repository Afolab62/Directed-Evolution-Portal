from __future__ import annotations

import logging
import sqlite3
from pathlib import Path
from typing import Any

from flask import current_app, has_app_context

from app.services.fingerprint import load_variants_from_db
from app.services.mutation_analysis import infer_uniprot_id, read_fasta

logger = logging.getLogger(__name__)

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
_INSTANCE_DB = _PROJECT_ROOT / "instance" / "app.db"

_WT_PROTEIN_CANDIDATES = (
    _PROJECT_ROOT / "data" / "wt_protein.fasta",
    _PROJECT_ROOT / "wt_protein.fasta",
    _PROJECT_ROOT / "O34996.fasta.txt",
)
_WT_PLASMID_CANDIDATES = (
    _PROJECT_ROOT / "data" / "pET-28a_BSU_DNA_Pol_I_WT.fa",
    _PROJECT_ROOT / "pET-28a_BSU_DNA_Pol_I_WT.fa",
    _PROJECT_ROOT / "data" / "wt_dna.fasta",
)


def _resolve_db_path(db_path: Path | None = None) -> Path:
    if db_path is not None:
        return Path(db_path)
    if has_app_context():
        return Path(current_app.config["DB_PATH"])
    return _INSTANCE_DB


def _connect(db_path: Path | None = None) -> sqlite3.Connection:
    conn = sqlite3.connect(_resolve_db_path(db_path))
    conn.row_factory = sqlite3.Row
    return conn


def _first_existing(paths: tuple[Path, ...]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def _ensure_wt_experiment_table(db_path: Path | None = None) -> None:
    conn = _connect(db_path)
    try:
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS wt_experiment (
                id INTEGER PRIMARY KEY,
                uniprot_id TEXT NOT NULL,
                wt_protein_seq TEXT NOT NULL,
                wt_plasmid_seq TEXT NOT NULL,
                updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
            """
        )
        conn.commit()
    finally:
        conn.close()


def get_all_variants_with_scores(db_path: Path | None = None) -> list[dict[str, Any]]:
    conn = _connect(db_path)
    try:
        rows = conn.execute(
            """
            SELECT
                Plasmid_Variant_Index,
                Parent_Plasmid_Variant,
                Directed_Evolution_Generation,
                Protein_Quantification_pg
            FROM activity_measurements
            WHERE Plasmid_Variant_Index IS NOT NULL
            ORDER BY Directed_Evolution_Generation, Plasmid_Variant_Index
            """
        ).fetchall()
    finally:
        conn.close()

    out: list[dict[str, Any]] = []
    for row in rows:
        try:
            out.append(
                {
                    "variant_id": int(row["Plasmid_Variant_Index"]),
                    "parent_variant_id": int(row["Parent_Plasmid_Variant"])
                    if row["Parent_Plasmid_Variant"] is not None
                    else -1,
                    "generation": int(row["Directed_Evolution_Generation"])
                    if row["Directed_Evolution_Generation"] is not None
                    else -1,
                    "activity_score": float(row["Protein_Quantification_pg"])
                    if row["Protein_Quantification_pg"] is not None
                    else None,
                }
            )
        except (TypeError, ValueError):
            continue
    return out


def save_active_experiment(
    uniprot_id: str,
    wt_protein_seq: str,
    wt_plasmid_seq: str,
    db_path: Path | None = None,
) -> None:
    _ensure_wt_experiment_table(db_path)
    conn = _connect(db_path)
    try:
        conn.execute(
            """
            INSERT INTO wt_experiment (id, uniprot_id, wt_protein_seq, wt_plasmid_seq, updated_at)
            VALUES (1, ?, ?, ?, CURRENT_TIMESTAMP)
            ON CONFLICT(id) DO UPDATE SET
                uniprot_id = excluded.uniprot_id,
                wt_protein_seq = excluded.wt_protein_seq,
                wt_plasmid_seq = excluded.wt_plasmid_seq,
                updated_at = CURRENT_TIMESTAMP
            """,
            (uniprot_id, wt_protein_seq, wt_plasmid_seq),
        )
        conn.commit()
        logger.info("Saved active experiment: uniprot_id=%s", uniprot_id)
    except Exception:
        conn.rollback()
        logger.exception("Could not save active experiment")
    finally:
        conn.close()


def get_staging_info(db_path: Path | None = None) -> dict[str, Any] | None:
    _ensure_wt_experiment_table(db_path)

    conn = _connect(db_path)
    try:
        row = conn.execute(
            """
            SELECT uniprot_id, wt_protein_seq, wt_plasmid_seq
            FROM wt_experiment
            WHERE id = 1
            """
        ).fetchone()
    finally:
        conn.close()

    if row and row["uniprot_id"] and row["wt_protein_seq"] and row["wt_plasmid_seq"]:
        logger.info("Loaded staging info from wt_experiment: %s", row["uniprot_id"])
        return {
            "uniprot_id": str(row["uniprot_id"]),
            "wt_protein_seq": str(row["wt_protein_seq"]),
            "wt_plasmid_seq": str(row["wt_plasmid_seq"]),
        }

    wt_dna_path = _first_existing(_WT_PLASMID_CANDIDATES)
    wt_protein_path = _first_existing(_WT_PROTEIN_CANDIDATES)
    if wt_dna_path is None or wt_protein_path is None:
        logger.warning("No staging row found and no WT FASTA files available")
        return None

    _, wt_plasmid_seq = read_fasta(wt_dna_path)
    wt_header, wt_protein_seq = read_fasta(wt_protein_path)
    uniprot_id = infer_uniprot_id(wt_header) or "O34996"

    if uniprot_id and wt_protein_seq and wt_plasmid_seq:
        save_active_experiment(
            uniprot_id=uniprot_id,
            wt_protein_seq=wt_protein_seq,
            wt_plasmid_seq=wt_plasmid_seq,
            db_path=db_path,
        )

    return {
        "wt_plasmid_seq": wt_plasmid_seq,
        "wt_protein_seq": wt_protein_seq,
        "uniprot_id": uniprot_id,
    }


__all__ = [
    "load_variants_from_db",
    "get_all_variants_with_scores",
    "get_staging_info",
    "save_active_experiment",
]