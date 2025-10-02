import argparse
from pathlib import Path

import pandas as pd


def normalize(series: pd.Series) -> pd.Series:
    series = series.fillna("").astype(str)
    return series.str.split(".").str[0]


def main() -> None:
    parser = argparse.ArgumentParser(description="Check ORF names from DeepBGC unique proteins against MIBiG protein IDs.")
    parser.add_argument(
        "--orf",
        default="checkmatch/final_output/unique_proteins_with_orf.tsv",
        help="Path to TSV created by the ORF annotator (columns: accession, orf_name)",
    )
    parser.add_argument(
        "--mibig",
        default="MIBiG_BGC/output/mibig_fasta_proteins_fixed.tsv",
        help="Path to MIBiG FASTA-derived TSV (columns include accession, protein_id)",
    )
    parser.add_argument(
        "--outdir",
        default="checkmatch/final_output",
        help="Directory to write the ORF match reports",
    )
    args = parser.parse_args()

    orf_path = Path(args.orf)
    mibig_path = Path(args.mibig)
    outdir = Path(args.outdir)

    if not orf_path.exists():
        raise SystemExit(f"ORF TSV not found: {orf_path}")
    if not mibig_path.exists():
        raise SystemExit(f"MIBiG TSV not found: {mibig_path}")

    outdir.mkdir(parents=True, exist_ok=True)

    orf_df = pd.read_csv(orf_path, sep="\t", dtype=str)
    mibig_df = pd.read_csv(mibig_path, sep="\t", dtype=str)

    required_orf_cols = {"accession", "orf_name"}
    required_mibig_cols = {"accession", "protein_id"}
    if not required_orf_cols.issubset(orf_df.columns):
        missing = ", ".join(sorted(required_orf_cols - set(orf_df.columns)))
        raise SystemExit(f"Missing column(s) in ORF TSV: {missing}")
    if not required_mibig_cols.issubset(mibig_df.columns):
        missing = ", ".join(sorted(required_mibig_cols - set(mibig_df.columns)))
        raise SystemExit(f"Missing column(s) in MIBiG TSV: {missing}")

    orf_df = orf_df.copy()
    mibig_df = mibig_df.copy()

    orf_df["acc_norm"] = normalize(orf_df["accession"])
    orf_df["orf_norm"] = normalize(orf_df["orf_name"])

    mibig_df["acc_norm"] = normalize(mibig_df["accession"])
    mibig_df["pid_norm"] = normalize(mibig_df["protein_id"])

    pid_lookup = (
        mibig_df.dropna(subset=["acc_norm", "pid_norm"])
        .groupby("acc_norm")["pid_norm"]
        .agg(lambda values: set(values))
    )

    def has_match(row: pd.Series) -> bool:
        candidates = pid_lookup.get(row["acc_norm"], set())
        return row["orf_norm"] in candidates

    orf_df["orf_in_mibig"] = orf_df.apply(has_match, axis=1)

    matches = orf_df.loc[orf_df["orf_in_mibig"]]
    missing = orf_df.loc[~orf_df["orf_in_mibig"]]

    matches_out = outdir / "orf_matches.tsv"
    missing_out = outdir / "orf_missing.tsv"
    full_out = outdir / "orf_match_results.tsv"

    matches.to_csv(matches_out, sep="\t", index=False)
    missing.to_csv(missing_out, sep="\t", index=False)
    orf_df.to_csv(full_out, sep="\t", index=False)

    print(f"Total ORF rows: {len(orf_df)}")
    print(f"Matches written to {matches_out} ({len(matches)})")
    print(f"Missing written to {missing_out} ({len(missing)})")
    print(f"Full results written to {full_out}")


if __name__ == "__main__":
    main()
