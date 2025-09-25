import argparse
import sys
from collections import Counter
from pathlib import Path
from typing import Dict

import pandas as pd

GROUPS = {
    "hydrophobic": set("AVILMFWY"),
    "polar": set("STNQ"),
    "positive": set("KRH"),
    "negative": set("DE"),
    "special": set("GPC"),
}


def aa_property_counts(seq: str, normalize: bool = False) -> Dict[str, float]:
    seq = (seq or "").strip().upper()
    counts = Counter(seq)
    total = sum(counts.values())
    group_counts: Dict[str, float] = {name: 0 for name in GROUPS}
    for group, aa_set in GROUPS.items():
        group_counts[group] = sum(counts.get(aa, 0) for aa in aa_set)
    if normalize and total > 0:
        group_counts = {group: value / total for group, value in group_counts.items()}
    return group_counts


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract physicochemical property group features from protein sequences"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to input TSV (with accession, protein_id, sequence columns)",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to output TSV for group features",
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
        props = aa_property_counts(row["sequence"], normalize=args.normalize)
        props["accession"] = row["accession"]
        props["protein_id"] = row["protein_id"]
        rows.append(props)

    columns = ["accession", "protein_id", *GROUPS.keys()]
    out_df = pd.DataFrame(rows, columns=columns)

    output_path = Path(args.output)
    if output_path.parent:
        output_path.parent.mkdir(parents=True, exist_ok=True)

    out_df.to_csv(output_path, sep="\t", index=False)

    print(
        f"Processed {len(df)} proteins. Wrote: {output_path} (normalize={args.normalize})"
    )


if __name__ == "__main__":
    main()
