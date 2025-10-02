import sys
import argparse
from pathlib import Path
import re
import pandas as pd


def normalize_accession_deepbgc(s: pd.Series) -> pd.Series:
    
    s = s.fillna("").astype(str)
    return s.str.split(".").str[0]


def normalize_accession_fasta(s: pd.Series) -> pd.Series:
    s = s.fillna("").astype(str)
    return s.str.split(".").str[0]


def normalize_protein_id(s: pd.Series) -> pd.Series:
    s = s.fillna("").astype(str)
    return s.str.split(".").str[0]


def acc_numeric(a: str) -> int:
    m = re.search(r"BGC(\d+)", a)
    return int(m.group(1)) if m else 10**12


def main():
    parser = argparse.ArgumentParser(description="Compare DeepBGC cleaned PFAM accessions with MIBiG FASTA accessions.")
    parser.add_argument(
        "--deepbgc",
        default="DeepBGC/Data/cleaned_deepbgc.tsv",
        help="Path to DeepBGC cleaned TSV (with 'accession' column)",
    )
    parser.add_argument(
        "--fasta",
        default="MIBiG_BGC/output/mibig_fasta_proteins_fixed.tsv",
        help="Path to MIBiG FASTA-derived TSV (with 'accession' column)",
    )
    parser.add_argument(
        "--outdir",
        default="checkmatch/final_output",
        help="Directory to write outputs",
    )
    args = parser.parse_args()

    deep_path = Path(args.deepbgc)
    fasta_path = Path(args.fasta)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not deep_path.exists():
        print(f"Error: DeepBGC TSV not found: {deep_path}")
        sys.exit(1)
    if not fasta_path.exists():
        print(f"Error: FASTA TSV not found: {fasta_path}")
        sys.exit(1)

    deep_df = pd.read_csv(deep_path, sep="\t", dtype=str)
    fasta_df = pd.read_csv(fasta_path, sep="\t", dtype=str)

    if "accession" not in deep_df.columns:
        print(f"Error: 'accession' column missing in {deep_path}")
        sys.exit(1)
    if "accession" not in fasta_df.columns:
        print(f"Error: 'accession' column missing in {fasta_path}")
        sys.exit(1)
    if "protein_id" not in deep_df.columns:
        print(f"Error: 'protein_id' column missing in {deep_path}")
        sys.exit(1)
    if "protein_id" not in fasta_df.columns:
        print(f"Error: 'protein_id' column missing in {fasta_path}")
        sys.exit(1)

    deep_df = deep_df.copy()
    fasta_df = fasta_df.copy()
    deep_df["acc_norm"] = normalize_accession_deepbgc(deep_df["accession"]) 
    deep_df["pid_norm"] = normalize_protein_id(deep_df["protein_id"]) 
    fasta_df["acc_norm"] = normalize_accession_fasta(fasta_df["accession"]) 
    fasta_df["pid_norm"] = normalize_protein_id(fasta_df["protein_id"]) 

    fasta_acc_set = set(fasta_df["acc_norm"].dropna())
    fasta_pid_set = set(fasta_df["pid_norm"].dropna())

    acc_unique_mask = ~deep_df["acc_norm"].isin(fasta_acc_set)
    pid_unique_mask = ~deep_df["pid_norm"].isin(fasta_pid_set)
    shared_acc_mask = ~acc_unique_mask

    unique_clusters = (
        deep_df.loc[acc_unique_mask, "acc_norm"]
        .dropna()
        .drop_duplicates()
        .sort_values(key=lambda s: s.map(acc_numeric))
    )
    unique_clusters_path = outdir / "unique_clusters.tsv"
    unique_clusters.to_frame(name="accession").to_csv(unique_clusters_path, sep="\t", index=False)

    unique_protein_map = (
        deep_df.loc[shared_acc_mask & pid_unique_mask, ["acc_norm", "pid_norm"]]
        .drop_duplicates()
        .rename(columns={"acc_norm": "accession", "pid_norm": "protein_id"})
    )
    unique_protein_map["source"] = "DeepBGC-protein-unique"
    unique_proteins_path = outdir / "unique_proteins.tsv"
    unique_protein_map.to_csv(unique_proteins_path, sep="\t", index=False)

    mibig_map = fasta_df[["acc_norm", "pid_norm"]].drop_duplicates().rename(columns={"acc_norm": "accession", "pid_norm": "protein_id"})
    mibig_map["source"] = "MIBiG"
    combined_map = pd.concat([mibig_map, unique_protein_map], ignore_index=True)
    combined_map = combined_map.assign(
        _acc_sort=combined_map["accession"].astype(str).map(acc_numeric)
    ).sort_values(by=['_acc_sort', 'accession']).drop(columns='_acc_sort').reset_index(drop=True)
    combined_map_path = outdir / "combined_unique_plus_mibig_mapping.tsv"
    combined_map.to_csv(combined_map_path, sep="\t", index=False)

    print(f"Unique DeepBGC clusters: {len(unique_clusters)} -> {unique_clusters_path}")
    print(f"Unique DeepBGC proteins (accession shared): {len(unique_protein_map)} -> {unique_proteins_path}")
    print(f"Combined mapping rows (MIBiG + DeepBGC protein-unique): {len(combined_map)} -> {combined_map_path}")


if __name__ == "__main__":
    main()
