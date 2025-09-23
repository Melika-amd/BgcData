import os
import json
import csv
import re
import argparse


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(description="Combine MIBiG JSON metadata into a TSV.")
    parser.add_argument(
        "--input",
        default=os.path.join(script_dir, "json_files", "mibig_json_4.0"),
        help="Path to folder containing MIBiG JSON files",
    )
    parser.add_argument(
        "--output",
        default=os.path.join(script_dir, "..", "output", "mibig_json_metadata.tsv"),
        help="Path to output TSV file",
    )
    args = parser.parse_args()

    folder_path = args.input
    output_file = args.output

    os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)

    headers = ["accession", "taxonomy", "compound", "bioactivities", "genome_accession", "location_start", "location_end"]

    # Collect all rows first so we can sort by accession numeric part
    rows = []
    for filename in os.listdir(folder_path):
        if filename.endswith(".json"):
            with open(os.path.join(folder_path, filename), "r") as file:
                data = json.load(file)

                accession = data.get("accession", "")
                taxonomy = data.get("taxonomy", {}).get("name", "")
                compound = data.get("compounds", [{}])[0].get("name", "") if data.get("compounds") else ""

                bioactivities = ""
                if data.get("compounds"):
                    bioactivity_list = []
                    for activity in data["compounds"][0].get("bioactivities", []):
                        name_field = activity.get("name")
                        if isinstance(name_field, dict):
                            bioactivity_list.append(name_field.get("activity", ""))
                        else:
                            bioactivity_list.append(str(name_field))
                    bioactivities = ", ".join(bioactivity_list)

                genome_accession = data.get("loci", [{}])[0].get("accession", "") if data.get("loci") else ""
                location_start = data.get("loci", [{}])[0].get("location", {}).get("from", "")
                location_end = data.get("loci", [{}])[0].get("location", {}).get("to", "")

                rows.append([
                    accession,
                    taxonomy,
                    compound,
                    bioactivities,
                    genome_accession,
                    location_start,
                    location_end
                ])

    # Sort by the numeric part of accession, e.g. BGC0000241 < BGC0000611
    def acc_num(row):
        acc = row[0]
        m = re.search(r"BGC(\d+)", acc)
        return int(m.group(1)) if m else float("inf")

    rows.sort(key=acc_num)

    with open(output_file, "w", newline='', encoding="utf-8") as tsv_file:
        writer = csv.writer(tsv_file, delimiter="\t")
        writer.writerow(headers)
        writer.writerows(rows)

    print(f"Done! Wrote: {output_file}")


if __name__ == "__main__":
    main()
