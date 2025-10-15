from collections import Counter
import csv


def main():
    input_path = "protein_id_types_with_accession.tsv"
    output_path = "protein_id_type_summary.tsv"

    counter = Counter()

    with open(input_path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            datatype = (row.get("datatype") or "").strip()
            if not datatype:
                datatype = "Unknown"
            counter[datatype] += 1

    total = sum(counter.values())
    if total == 0:
        raise RuntimeError(f"No rows found in {input_path}")

    with open(output_path, "w", encoding="utf-8", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t")
        writer.writerow(["datatype", "count", "percentage"])
        for datatype, count in counter.most_common():
            percentage = (count / total) * 100
            writer.writerow([datatype, count, f"{percentage:.2f}"])
        writer.writerow(["TOTAL", total, "100.00"])

    non_ncbi_count = counter.get("Not in NCBI", 0)
    non_ncbi_pct = (non_ncbi_count / total) * 100
    print(f" Wrote summary to {output_path}")
    print(f" 'Not in NCBI' entries: {non_ncbi_count} ({non_ncbi_pct:.2f}%)")


if __name__ == "__main__":
    main()
