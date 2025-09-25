import argparse
import os
import sys
from collections import Counter
from pathlib import Path
from typing import Dict

import pandas as pd

AA20 = list("ACDEFGHIKLMNPQRSTVWY")


def aa_composition(seq: str, normalize: bool = False) -> Dict[str, float]:
    seq = (seq or "").strip().upper()
    counts = Counter(seq)
    total = sum(counts[aa] for aa in AA20)
    comp = {aa: counts.get(aa, 0) for aa in AA20}
    if normalize and total > 0:
        comp = {aa: val / total for aa, val in comp.items()}
    return comp


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract amino acid composition from protein sequences"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to input TSV (with accession, protein_id, sequence columns)",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to output TSV for amino acid composition",
    )
    parser.add_argument(
        "--normalize",
        action="store_true",
        help="Normalize counts to frequencies (0–1)",
    )
    args = parser.parse_args()

    try:
        df = pd.read_csv(
            args.input,
            sep="\t",
            dtype=str,
            usecols=["accession", "protein_id", "sequence"],
        )
    except FileNotFoundError:
        print(f"Input file not found: {args.input}")
        sys.exit(1)

    for col in ("accession", "protein_id", "sequence"):
        if col not in df.columns:
            print(f"Missing required column in input: {col}")
            sys.exit(1)

    rows = []
    for _, row in df.iterrows():
        comp = aa_composition(row["sequence"], normalize=args.normalize)
        comp["accession"] = row["accession"]
        comp["protein_id"] = row["protein_id"]
        rows.append(comp)

    columns = ["accession", "protein_id", *AA20]
    out_df = pd.DataFrame(rows, columns=columns)

    output_path = Path(args.output)
    if output_path.parent:
        os.makedirs(output_path.parent, exist_ok=True)

    out_df.to_csv(output_path, sep="\t", index=False)

    print(
        f"Processed {len(df)} proteins. Wrote: {output_path} (normalize={args.normalize})"
    )


if __name__ == "__main__":
    main()
