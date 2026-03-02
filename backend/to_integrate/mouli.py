# %%
import pandas as pd
from datetime import datetime

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from sqlalchemy.orm import declarative_base, relationship
from sqlalchemy import Column, Integer, Float, String, Boolean, ForeignKey

# %%
# SQLite database file (local)
DB_FILE = "samples.db"
engine = create_engine(f"sqlite:///{DB_FILE}", echo=True)

Session = sessionmaker(bind=engine)
session = Session()

# %%
Base = declarative_base()

class Sample(Base):
    __tablename__ = "samples"

    id = Column(Integer, primary_key=True, autoincrement=True)
    plasmid_variant_index = Column(Float, nullable=False)
    parent_plasmid_variant = Column(Float, nullable=False)
    directed_evolution_generation = Column(Float, nullable=False)
    assembled_dna_sequence = Column(String(2000), nullable=False)
    dna_quantification_fg = Column(Float, nullable=False)
    protein_quantification_pg = Column(Float, nullable=False)
    is_control = Column(Boolean, nullable=False)

    extra_data = relationship("SampleMetadata", back_populates="sample", cascade="all, delete-orphan")


class SampleMetadata(Base):
    __tablename__ = "sample_metadata"

    id = Column(Integer, primary_key=True, autoincrement=True)
    sample_id = Column(Integer, ForeignKey("samples.id"), nullable=False)
    key = Column(String(255), nullable=False)
    value = Column(String(2000))

    sample = relationship("Sample", back_populates="extra_data")

# %%
Base.metadata.create_all(engine)

# %%
essential_fields = {
    "plasmid_variant_index": float,
    "parent_plasmid_variant": float,
    "directed_evolution_generation": float,

    "assembled_dna_sequence": str,
    "dna_quantification_fg": float,
    "protein_quantification_pg": float,
    "is_control": bool
}

col_synonyms = {
    "plasmid_variant_index": ["variant_index", "plasmid_id"],
    "parent_plasmid_variant": ["parent_variant", "parent_id"],
    "directed_evolution_generation": ["generation", "evolution_generation"],
    "assembled_dna_sequence": ["dna_sequence", "sequence", "assembled_sequence"],
    "dna_quantification_fg": ["dna_concentration_fg", "dna_qty_fg"],
    "protein_quantification_pg": ["protein_concentration_pg", "protein_qty_pg"],
    "is_control": ["control", "is_control", "control_sample"]
}

# %%
def clean_cols(col: str) -> str:
    return col.strip().lower().replace(" ", "_").replace("-", "_")

def build_synonym_map(col_synonyms): #Reverse synonym lookup
    synonym_map = {}
    for synonym, variants in col_synonyms.items():
        for v in variants:
            synonym_map[clean_cols(v)] = synonym
    return synonym_map

def validate_mapping(mapping): #Prevents uplicate assignments. Remove if confirmation can be require in fronten. 
    reverse = {}
    for raw, field in mapping.items():
        if field in reverse:
            raise ValueError(
                f"Multiple columns mapped to '{field}': "
                f"{reverse[field]} and {raw}"
            )
        reverse[field] = raw

# %%
def confirm_mapping_bulk(mapping):

    while True:
        print("\nProposed Column Mapping:")
        for raw, canon in mapping.items():
            print(f"{raw} → {canon}")

        response = input("\nPress Enter to accept, or type 'edit' to modify: ").strip().lower()
        if response == "":
            break  # User accepted mapping

        if response == "edit":
            edits = input(
                "Enter edits as raw:target,raw2:target2,... : "
            ).strip()
            for pair in edits.split(","):
                if ":" not in pair:
                    print(f"Skipping invalid entry '{pair}'")
                    continue
                raw_col, target_field = pair.split(":", 1)
                raw_col = raw_col.strip()
                target_field = target_field.strip()
                if raw_col not in mapping:
                    print(f"Column '{raw_col}' not in proposed mapping. Skipping.")
                    continue
                mapping[raw_col] = target_field
            # After edits, loop prints updated mapping automatically
        else:
            print("Invalid input. Press Enter to accept or type 'edit' to modify.")

    return mapping

# %%
class ColumnMapper:
    def __init__(self, essential_fields, col_synonyms):
        self.essential_fields = list(essential_fields.keys())
        self.synonym_map = build_synonym_map(col_synonyms)

    def auto_map_by_synonym(self, columns):
        mapping = {}
        used_cols = set()

        for col in columns:
            if col in self.synonym_map:
                official = self.synonym_map[col]
                if official not in mapping.values():  # avoid duplicates
                    mapping[col] = official
                    used_cols.add(col)

        return mapping, used_cols

    def left_to_right_assign(self, columns, used_cols, existing_mapping): #If NOT already mapped

        remaining_cols = [c for c in columns if c not in used_cols]
        already_assigned_fields = set(existing_mapping.values())
        remaining_fields = [f for f in self.essential_fields if f not in already_assigned_fields]

        for col, field in zip(remaining_cols, remaining_fields):
            existing_mapping[col] = field

        return existing_mapping

    def generate_mapping(self, df_columns):
        # Track original ↔ cleaned names
        original_to_clean = {c: clean_cols(c) for c in df_columns}
        clean_to_original = {v: k for k, v in original_to_clean.items()}

        cleaned_cols = list(clean_to_original.keys())

        #Synonym mapping
        mapping, used = self.auto_map_by_synonym(cleaned_cols)

        #Left-to-right
        mapping = self.left_to_right_assign(cleaned_cols, used, mapping)

        #Validate before user sees it
        validate_mapping(mapping)

        mapped_fields = set(mapping.values())
        missing_fields = [f for f in self.essential_fields if f not in mapped_fields]

        if missing_fields:
            print("\nMissing essential fields (not found in file):", missing_fields)


        #User confirmation - do in front en later?
        mapping = confirm_mapping_bulk(mapping)

        # Convert cleaned names back to original DataFrame column names
        final_mapping = {clean_to_original[k]: v for k, v in mapping.items()}
        return final_mapping

# %%
def validate_row(row):
    errors = []

    # Missing essential fields
    for field in essential_fields:
        if pd.isna(row.get(field)):
            errors.append(f"Missing value for {field}")

    # Biological / logical rules
    if row.get("directed_evolution_generation", 0) < 0:
        errors.append("Generation cannot be negative")
    if row.get("dna_quantification_fg", 0) < 0:
        errors.append("DNA quantification cannot be negative")
    if row.get("protein_quantification_pg", 0) < 0:
        errors.append("Protein quantification cannot be negative")
    if row.get("is_control") not in [True, False]:
        errors.append("is_control must be boolean")

    return errors


# %%
def insert_sql(df_valid):
    metadata_columns = [c for c in df_valid.columns if c not in essential_fields]

    for _, row in df_valid.iterrows():
        sample = Sample(
            plasmid_variant_index=row["plasmid_variant_index"],
            parent_plasmid_variant=row["parent_plasmid_variant"],
            directed_evolution_generation=row["directed_evolution_generation"],
            assembled_dna_sequence=row["assembled_dna_sequence"],
            dna_quantification_fg=row["dna_quantification_fg"],
            protein_quantification_pg=row["protein_quantification_pg"],
            is_control=row["is_control"]
        )

        for col in metadata_columns:
            val = row[col]
            if pd.notna(val):
                sample.metadata.append(SampleMetadata(key=col, value=str(val)))

        session.add(sample)
    session.commit()


# %%
class FileLoader:
    def load(self, filepath) -> pd.DataFrame:
        if filepath.endswith(".tsv"):
            return pd.read_csv(filepath, sep="\t")
        elif filepath.endswith(".json"):
            return pd.read_json(filepath)
        else:
            raise ValueError("Unsupported format")


# %%
def coerce_types(df, essential_fields):
    df = df.copy()
    for col, dtype in essential_fields.items():
        if col not in df.columns:
            continue

        if dtype in [float, int]:
            df[col] = pd.to_numeric(df[col], errors="coerce")

        elif dtype == bool:
            df[col] = df[col].map(
                {1: True, 0: False, "1": True, "0": False,
                 "true": True, "false": False,
                 "True": True, "False": False}
            )

        elif dtype == str:
            df[col] = df[col].astype(str)

    return df

# %%
def validate_row(row):
    errors = []

    # Missing essential values
    for field in essential_fields:
        if pd.isna(row.get(field)):
            errors.append(f"Missing value for {field}")

    # Logical rules
    if pd.notna(row.get("directed_evolution_generation")) and row["directed_evolution_generation"] < 0:
        errors.append("Generation cannot be negative")

    if pd.notna(row.get("dna_quantification_fg")) and row["dna_quantification_fg"] < 0:
        errors.append("DNA quantification cannot be negative")

    if pd.notna(row.get("protein_quantification_pg")) and row["protein_quantification_pg"] < 0:
        errors.append("Protein quantification cannot be negative")

    if row.get("is_control") not in [True, False]:
        errors.append("is_control must be boolean")

    seq = row.get("assembled_dna_sequence")
    if isinstance(seq, str) and not set(seq.upper()).issubset(set("ATCG")):
        errors.append("DNA sequence contains invalid characters")

    return errors

# %%
def insert_sql(valid_data):
    session = Session()

    metadata_columns = [c for c in valid_data.columns if c not in essential_fields]

    for _, row in valid_data.iterrows():
        sample = Sample(
            plasmid_variant_index=row["plasmid_variant_index"],
            parent_plasmid_variant=row["parent_plasmid_variant"],
            directed_evolution_generation=row["directed_evolution_generation"],
            assembled_dna_sequence=row["assembled_dna_sequence"],
            dna_quantification_fg=row["dna_quantification_fg"],
            protein_quantification_pg=row["protein_quantification_pg"],
            is_control=row["is_control"]
        )

        # Store non-essential columns in metadata table
        for col in metadata_columns:
            val = row[col]
            if pd.notna(val):
                meta_entry = SampleMetadata(key=col, value=str(val))
                sample.metadata.append(meta_entry)

        session.add(sample)

    session.commit()
    session.close()

# %% [markdown]
# ### Pipeline

# %%
file_path = "c://Users//Leora//OneDrive - Queen Mary, University of London//Group_Project//Example_Data//DE_BSU_Pol_Batch_1.tsv"

loader = FileLoader()
df = loader.load(file_path)

mapper = ColumnMapper(essential_fields, col_synonyms)
column_mapping = mapper.generate_mapping(df.columns)
df = df.rename(columns=column_mapping)


# Ensure all essential columns exist
missing_fields = [f for f in essential_fields if f not in df.columns]
if missing_fields:
    raise ValueError(
        f"Missing essential fields in file: {missing_fields}. "
        "Can not continue until these columns are present."
    )

df = coerce_types(df, essential_fields)

valid_rows, rejected_rows = [], []

for idx, row in df.iterrows():
    errs = validate_row(row)
    if errs:
        row["qc_error_reason"] = "; ".join(errs)
        row["qc_row_number"] = idx + 1  # Human-readable row number
        rejected_rows.append(row)
    else:
        valid_rows.append(row)

df_valid = pd.DataFrame(valid_rows)
df_rejected = pd.DataFrame(rejected_rows)


print(f"\nQC Report:")
print(f"Valid rows: {len(df_valid)}")
print(f"Rejected rows: {len(df_rejected)}")

if len(df_rejected) > 0:
    print("\n Rejected Row Summary:")
    for _, r in df_rejected.iterrows():
        print(f"Row {r['qc_row_number']}: {r['qc_error_reason']}")

    raise ValueError(" QC Failed — Fix and reupload file.")

input("\n QC complete. Press Enter to continue.")


insert_sql(df_valid)



