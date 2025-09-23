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
    core = s.str.split("|").str[0]
    return core.str.split(".").str[0]


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
        default="MIBiG_BGC/output/mibig_fasta_proteins.tsv",
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

    deep_ids = set(normalize_accession_deepbgc(deep_df["accession"]))
    fasta_ids = set(normalize_accession_fasta(fasta_df["accession"]))

    unique_deep = sorted(deep_ids - fasta_ids, key=acc_numeric)

    
    unique_path = outdir / "unique_deepbgc_accessions.tsv"
    pd.DataFrame({"accession": unique_deep}).to_csv(unique_path, sep="\t", index=False)

    
    combined_records = (
        [(a, "MIBiG") for a in sorted(fasta_ids, key=acc_numeric)] +
        [(a, "DeepBGC-unique") for a in unique_deep]
    )
    combined_df = pd.DataFrame(combined_records, columns=["accession", "source"])
    combined_path = outdir / "combined_unique_plus_mibig.tsv"
    combined_df.to_csv(combined_path, sep="\t", index=False)

    print(f"Unique DeepBGC accessions: {len(unique_deep)} -> {unique_path}")
    print(f"Combined (MIBiG + DeepBGC-unique): {len(combined_df)} -> {combined_path}")


if __name__ == "__main__":
    main()
