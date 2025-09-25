import argparse
import csv
from collections import Counter
from pathlib import Path
from typing import Iterable, Tuple, List

from Bio import SeqIO


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract k-mer frequencies from a FASTA file and write them to a TSV table."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to input FASTA file",
    )
    parser.add_argument(
        "--k",
        type=int,
        required=True,
        help="K-mer length (positive integer)",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to output TSV file",
    )
    return parser.parse_args()


def get_kmers(sequence: str, k: int) -> Iterable[str]:
    sequence = sequence.upper()
    for i in range(len(sequence) - k + 1):
        yield sequence[i : i + k]


def collect_counts(fasta_path: Path, k: int) -> Tuple[List[str], List[Tuple[str, Counter]]]:
    all_kmers = set()
    per_sequence: List[Tuple[str, Counter]] = []

    for record in SeqIO.parse(str(fasta_path), "fasta"):
        seq = str(record.seq)
        if len(seq) < k:
            counts = Counter()
        else:
            counts = Counter(get_kmers(seq, k))
            all_kmers.update(counts.keys())
        per_sequence.append((record.id, counts))

    return sorted(all_kmers), per_sequence


def write_counts(output_path: Path, kmers: List[str], per_sequence: List[Tuple[str, Counter]]) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sequence_id", *kmers])

        for sequence_id, counts in per_sequence:
            row = [sequence_id]
            row.extend(counts.get(kmer, 0) for kmer in kmers)
            writer.writerow(row)


def main() -> None:
    args = parse_args()

    k = args.k
    if k <= 0:
        raise ValueError("--k must be a positive integer")

    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.exists():
        raise FileNotFoundError(f"Input FASTA not found: {input_path}")

    kmers, per_sequence = collect_counts(input_path, k)

    if not kmers:
        print("Warning: No k-mers were generated (sequences shorter than k or empty file).")

    write_counts(output_path, kmers, per_sequence)

    print(f"Done! Wrote k-mer counts to {output_path}")


if __name__ == "__main__":
    main()
