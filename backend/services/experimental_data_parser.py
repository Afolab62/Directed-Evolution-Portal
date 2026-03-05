"""
Service for parsing experimental data files (TSV/JSON).

Pipeline (mirrors mouli.py):
  1. Load file into DataFrame
  2. Map columns: synonym lookup then left-to-right positional fallback
  3. Coerce types (bool via explicit map, numeric via pd.to_numeric)
  4. Validate each row (soft QC — rejected rows are reported, not discarded)
  5. Split into valid_df / control_df / rejected_df
"""

import pandas as pd
import io
from typing import Dict, List, Any


# ── Canonical field names and DB column types ────────────────────────────────

# Required fields (upload is blocked only if these cannot be found at all)
REQUIRED_FIELDS = {
    "plasmid_variant_index": float,
    "generation": int,
    "assembled_dna_sequence": str,
    "dna_yield": float,
    "protein_yield": float,
    "is_control": bool,
}

# Optional fields (present in most files; absent for gen-0 variants)
OPTIONAL_FIELDS = {
    "parent_plasmid_variant": float,
}

ESSENTIAL_FIELDS = {**REQUIRED_FIELDS, **OPTIONAL_FIELDS}

# ── Synonyms for every canonical name ────────────────────────────────────────

COLUMN_SYNONYMS: Dict[str, List[str]] = {
    "plasmid_variant_index": [
        "plasmid_variant_index", "variant_index", "plasmid_id",
        "variant_id", "index",
    ],
    "parent_plasmid_variant": [
        "parent_plasmid_variant", "parent_variant", "parent_id",
        "parent", "parent_index",
    ],
    "generation": [
        "generation", "directed_evolution_generation",
        "evolution_generation", "gen",
    ],
    "assembled_dna_sequence": [
        "assembled_dna_sequence", "dna_sequence", "sequence",
        "assembled_sequence", "plasmid_sequence",
    ],
    "dna_yield": [
        "dna_yield", "dna_quantification_fg", "dna_qty_fg",
        "dna_concentration_fg", "dna_quantification",
    ],
    "protein_yield": [
        "protein_yield", "protein_quantification_pg", "protein_qty_pg",
        "protein_concentration_pg", "protein_quantification",
    ],
    "is_control": [
        "is_control", "control", "control_sample", "iscontrol",
    ],
}


def _clean(col: str) -> str:
    return col.strip().lower().replace(" ", "_").replace("-", "_")


def _build_synonym_map() -> Dict[str, str]:
    m: Dict[str, str] = {}
    for canonical, aliases in COLUMN_SYNONYMS.items():
        for alias in aliases:
            m[_clean(alias)] = canonical
    return m


_SYNONYM_MAP = _build_synonym_map()

# Bool coercion map 
_BOOL_MAP = {
    1: True, 0: False,
    "1": True, "0": False,
    "true": True, "false": False,
    "True": True, "False": False,
    "TRUE": True, "FALSE": False,
    "yes": True, "no": False,
    "Yes": True, "No": False,
    "YES": True, "NO": False,
}


# ── Row-level QC ─────────────────────────────

def _validate_row(row: Dict[str, Any]) -> List[str]:
    errors: List[str] = []

    # Required numeric fields must not be NaN
    for field in ("plasmid_variant_index", "generation", "dna_yield", "protein_yield"):
        if pd.isna(row.get(field)):
            errors.append(f"Missing value for '{field}'")

    if pd.notna(row.get("generation")) and row["generation"] < 0:
        errors.append("Generation cannot be negative")

    if pd.notna(row.get("dna_yield")) and row["dna_yield"] < 0:
        errors.append("dna_yield cannot be negative")

    if pd.notna(row.get("protein_yield")) and row["protein_yield"] < 0:
        errors.append("protein_yield cannot be negative")

    if row.get("is_control") not in (True, False):
        errors.append("is_control must be boolean (TRUE/FALSE/1/0)")

    seq = row.get("assembled_dna_sequence")
    if isinstance(seq, str) and seq and not set(seq.upper()).issubset(set("ATCGNRYZ")):
        errors.append("DNA sequence contains invalid characters (only A/T/C/G allowed)")

    return errors


# ── Main parser class ─────────────────────────────────────────────────────────

class ExperimentalDataParser:

    def __init__(self):
        self._all_field_names = list(ESSENTIAL_FIELDS.keys())
        self._required_field_names = list(REQUIRED_FIELDS.keys())

    def preview_mapping(self, file_content: str, file_format: str) -> dict:
        """
        Parse only the header of *file_content* and return the auto-detected
        column mapping together with the list of all raw column names.
        Does NOT perform type coercion, QC, or any database work.

        Returns:
            {
              "raw_columns": [<original column names> ...],
              "column_mapping": {<raw_col>: <canonical_field> ...},
              "metadata_columns": [<unmapped raw column names> ...],
              "missing_required": [<canonical fields not found> ...],
              "canonical_fields": [<all canonical field names> ...],
            }
        """
        df = self._parse(file_content, file_format)
        raw_columns = df.columns.tolist()
        column_mapping, missing_required = self._map_columns(raw_columns)
        assigned_raws = set(column_mapping.keys())
        metadata_columns = [c for c in raw_columns if c not in assigned_raws]
        return {
            "raw_columns": raw_columns,
            "column_mapping": column_mapping,
            "metadata_columns": metadata_columns,
            "missing_required": missing_required,
            "canonical_fields": self._all_field_names,
        }

    def process_file(self, file_content: str, file_format: str,
                     column_mapping_override: dict | None = None):
        """
        Parse a TSV or JSON upload and return:
            valid_df, control_df, rejected_df, summary

        Every parseable row is returned — rows that fail QC land in
        rejected_df but the upload is never blocked on their account.
        The upload is only blocked if the file itself cannot be read or
        if a required column cannot be found at all.

        If *column_mapping_override* is provided (a {raw_col: canonical_field}
        dict coming from the frontend mapping-confirmation step), it is used
        instead of auto-detection.  Auto-detection is still run to verify that
        all required fields are covered.
        """
        df = self._parse(file_content, file_format)

        # ── 5_change Detect duplicate rows across all columns ─────────────
        duplicates = df[df.duplicated(keep=False)]
        if not duplicates.empty:
            print("\nDuplicate rows detected:")
            for idx in duplicates.index:
                print(f"Row {idx + 1} is duplicated")
            raise ValueError(
                "Duplicate rows detected in input file. "
                "Remove duplicates before ingestion."
            )

        # ── Column mapping ────────────────────────────────────────────────
        if column_mapping_override:
            column_mapping = {k: v for k, v in column_mapping_override.items()
                              if k in df.columns and v}  # strip blanks / phantom cols
            assigned_canonicals = set(column_mapping.values())
            missing_required = [f for f in self._required_field_names
                                 if f not in assigned_canonicals]
        else:
            column_mapping, missing_required = self._map_columns(df.columns.tolist())

        if missing_required:
            raise ValueError(
                f"Could not find required columns: {', '.join(missing_required)}. "
                "Please check your file contains: "
                + ", ".join(self._required_field_names)
            )

        df = df.rename(columns=column_mapping)
        metadata_columns = [c for c in df.columns if c not in self._all_field_names]

        # ── Type coercion ─────────────────────────────────────────────────
        df = self._coerce(df)

        # ── Row-level QC + split ──────────────────────────────────────────
        valid_rows: List[Dict] = []
        control_rows: List[Dict] = []
        rejected_rows: List[Dict] = []

        for idx, row in df.iterrows():
            d = row.to_dict()
            errors = _validate_row(d)
            if errors:
                d["qc_error_reason"] = "; ".join(errors)
                d["qc_row_number"] = int(idx) + 1
                rejected_rows.append(d)
            elif d.get("is_control") is True:
                control_rows.append(d)
            else:
                valid_rows.append(d)

        valid_df = pd.DataFrame(valid_rows) if valid_rows else pd.DataFrame()
        control_df = pd.DataFrame(control_rows) if control_rows else pd.DataFrame()
        rejected_df = pd.DataFrame(rejected_rows) if rejected_rows else pd.DataFrame()

        summary = {
            "total_rows": len(df),
            "valid_rows": len(valid_df),
            "control_rows": len(control_df),
            "rejected_rows": len(rejected_df),
            "column_mapping": column_mapping,
            "metadata_columns": metadata_columns,
            "rejected_details": [
                {"qc_row_number": r["qc_row_number"], "qc_error_reason": r["qc_error_reason"]}
                for r in rejected_rows
            ],
        }
        return valid_df, control_df, rejected_df, summary

    # ── Internal helpers ──────────────────────────────────────────────────────

    def _parse(self, content: str, fmt: str) -> pd.DataFrame:
        try:
            if fmt.lower() in ("tsv", "txt"):
                df = pd.read_csv(io.StringIO(content), sep="\t")
            elif fmt.lower() == "json":
                df = pd.read_json(io.StringIO(content))
            else:
                raise ValueError(f"Unsupported format: {fmt}")
            if df.empty:
                raise ValueError("File contains no data rows")
            return df
        except Exception as exc:
            raise ValueError(f"Failed to parse file: {exc}") from exc

    def _map_columns(self, df_columns: List[str]):
        """
        Pass 1: synonym lookup (cleaned names).
        Pass 2: left-to-right positional assignment for still-unmapped columns.

        Returns (mapping_dict, missing_required_fields).
        """
        original_to_clean = {c: _clean(c) for c in df_columns}

        mapping: Dict[str, str] = {}
        used_originals: set = set()
        assigned_canonicals: set = set()

        # Pass 1 — synonym lookup
        for orig, cleaned in original_to_clean.items():
            canonical = _SYNONYM_MAP.get(cleaned)
            if canonical and canonical not in assigned_canonicals:
                mapping[orig] = canonical
                used_originals.add(orig)
                assigned_canonicals.add(canonical)

        # Pass 2 — positional fallback ONLY for unmatched *required* fields.
        # Optional fields (e.g. parent_plasmid_variant) are intentionally excluded:
        # if they can't be found by name they should stay absent, not silently
        # consume the first unrecognised / extra-metadata column.
        remaining_originals = [c for c in df_columns if c not in used_originals]
        remaining_canonicals = [
            f for f in self._required_field_names if f not in assigned_canonicals
        ]
        for orig, canonical in zip(remaining_originals, remaining_canonicals):
            mapping[orig] = canonical
            assigned_canonicals.add(canonical)

        missing_required = [f for f in self._required_field_names if f not in assigned_canonicals]
        return mapping, missing_required

    def _coerce(self, df: pd.DataFrame) -> pd.DataFrame:
        """Coerce each essential column to its target type."""
        df = df.copy()
        for col, dtype in ESSENTIAL_FIELDS.items():
            if col not in df.columns:
                continue
            if dtype in (float, int):
                df[col] = pd.to_numeric(df[col], errors="coerce")
                if dtype == int:
                    df[col] = df[col].astype("Int64")
            elif dtype == bool:
                # Use explicit map (mirrors mouli.py) — anything not in the
                # map becomes NaN, which is caught by _validate_row
                df[col] = df[col].map(_BOOL_MAP)
            elif dtype == str:
                # Preserve NaN — don't stringify it
                df[col] = df[col].where(df[col].isna(), df[col].astype(str))
        return df


# Singleton used by the route layer
parser = ExperimentalDataParser()
