import pandas as pd
import os

print(" Current Working Directory:", os.getcwd())

# Resolve paths relative to this script, so CWD doesn't matter
script_dir = os.path.dirname(os.path.abspath(__file__))
deepbgc_dir = os.path.dirname(script_dir)
data_dir = os.path.join(deepbgc_dir, "Data")
input_path = os.path.join(data_dir, "MIBiG.pfam.tsv")
output_path = os.path.join(data_dir, "cleaned_deepbgc.tsv")

# Load only the MIBiG PFAM data
try:
    mibig_df = pd.read_csv(input_path, sep="\t")
except FileNotFoundError as e:
    print(f" File not found: {e}")
    exit(1)

# Select and rename columns
required_cols = [
    'sequence_id', 'protein_id', 'gene_start', 'gene_end', 'gene_strand', 'pfam_id', 'in_cluster'
]
missing = [c for c in required_cols if c not in mibig_df.columns]
if missing:
    print(f" Missing required columns in MIBiG.pfam.tsv: {missing}")
    exit(1)

# Keep only relevant columns and rename sequence_id -> accession
clean_df = (
    mibig_df[required_cols]
    .rename(columns={'sequence_id': 'accession'})
)

# Remove version suffix like ".1", ".2" from accession (e.g., BGC0000001.1 -> BGC0000001)
clean_df['accession'] = clean_df['accession'].astype(str).str.replace(r'\.\d+$', '', regex=True)

# Basic cleaning: drop rows with missing key fields, drop duplicates
clean_df = clean_df.dropna(subset=['accession', 'protein_id', 'pfam_id'])
clean_df = clean_df.drop_duplicates()

# Write cleaned dataset
clean_df.to_csv(output_path, sep="\t", index=False)
print(f"\n Saved cleaned dataset to: {output_path} ({len(clean_df):,} rows)")
