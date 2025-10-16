
"""
Produce a cleaned TSV from `unique_proteins_with_genbank.tsv` that keeps only
the columns `accession`, `protein_id`, and `datatype`, while normalising the
values for entries that were previously marked as "Not in NCBI".

Rules applied per row:
- Existing GenBank rows remain unchanged.
- "Not in NCBI" rows with a `new_genbank_id` use that value as the `protein_id`
  and their `datatype` is set to "GenBank".
- "Not in NCBI" rows without a resolved ID keep their original `protein_id`
  and have `datatype` set to "Unknown".
"""
from __future__ import annotations

import csv
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
SOURCE_FILE = SCRIPT_DIR / "unique_proteins_with_genbank.tsv"
OUTPUT_FILE = SCRIPT_DIR / "unique_proteins_normalized.tsv"


def resolve_values(row: dict[str, str]) -> tuple[str, str]:
    datatype = (row.get("datatype") or "").strip()
    protein_id = (row.get("protein_id") or "").strip()
    new_id = (row.get("new_genbank_id") or "").strip()

    lower_type = datatype.lower()
    if "genbank" in lower_type:
        return protein_id, "GenBank"
    if "not in ncbi" in lower_type:
        if new_id:
            return new_id, "GenBank"
        return protein_id, "Unknown"
    return protein_id or "Unknown", datatype or "Unknown"


def main() -> int:
    if not SOURCE_FILE.exists():
        raise FileNotFoundError(f"Source TSV not found: {SOURCE_FILE}")

    with SOURCE_FILE.open("r", newline="", encoding="utf-8") as src, OUTPUT_FILE.open(
        "w", newline="", encoding="utf-8"
    ) as dst:
        reader = csv.DictReader(src, delimiter="\t")
        required = {"accession", "protein_id", "datatype", "new_genbank_id"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Input file missing columns: {', '.join(sorted(missing))}")

        fieldnames = ["accession", "protein_id", "datatype"]
        writer = csv.DictWriter(dst, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for row in reader:
            protein_id, datatype = resolve_values(row)
            writer.writerow(
                {
                    "accession": (row.get("accession") or "").strip(),
                    "protein_id": protein_id,
                    "datatype": datatype,
                }
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

