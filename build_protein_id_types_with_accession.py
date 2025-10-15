def main():
    proteins_path = "MIBiG_BGC/output/mibig_fasta_proteins_fixed.tsv"
    types_path = "protein_id_types.tsv"
    output_path = "protein_id_types_with_accession.tsv"

    id_type_map = {}
    with open(types_path, "r", encoding="utf-8") as handle:
        next(handle)  
        for row in handle:
            row = row.strip()
            if not row:
                continue
            protein_id, id_type = row.split("\t", 1)
            id_type_map[protein_id] = id_type

    seen = set()
    count = 0

    with open(proteins_path, "r", encoding="utf-8") as proteins, open(
        output_path, "w", encoding="utf-8"
    ) as out:
        header = proteins.readline().strip().split("\t")
        try:
            accession_idx = header.index("accession")
            protein_idx = header.index("protein_id")
        except ValueError as exc:
            raise RuntimeError("Expected 'accession' and 'protein_id' columns") from exc

        out.write("accession\tprotein_id\tdatatype\n")

        for line in proteins:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(accession_idx, protein_idx):
                continue

            accession = parts[accession_idx]
            protein_id = parts[protein_idx]

            if not protein_id:
                continue

            key = (accession, protein_id)
            if key in seen:
                continue
            seen.add(key)

            datatype = id_type_map.get(protein_id, "Unknown")
            out.write(f"{accession}\t{protein_id}\t{datatype}\n")
            count += 1

    print(f" Wrote {count} rows to {output_path}")


if __name__ == "__main__":
    main()
