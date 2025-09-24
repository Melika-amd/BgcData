from Bio import SeqIO
import os
import csv
import argparse


def main():
    parser = argparse.ArgumentParser(description="Extract protein sequences from MIBiG FASTA to TSV.")
    parser.add_argument(
        "--input",
        default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "FASTAFiles", "mibig_prot_seqs_4.0.fasta"),
        help="Path to input FASTA file",
    )
    parser.add_argument(
        "--output",
        default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "output", "mibig_fasta_proteins_fixed.tsv"),
        help="Path to output TSV file",
    )
    args = parser.parse_args()

    input_fasta = args.input
    output_tsv = args.output

 
    os.makedirs(os.path.dirname(os.path.abspath(output_tsv)), exist_ok=True)

    with open(output_tsv, "w", newline="", encoding="utf-8") as tsv_file:
        writer = csv.writer(tsv_file, delimiter="\t")
      
        writer.writerow([
            "accession",   # part[0]
            "protein_id",  # part[4]
            "contig",      # part[1]
            "location",    # part[2]
            "strand",      # part[3]
            "description", # part[5]
            "sequence",
        ])

        for record in SeqIO.parse(input_fasta, "fasta"):
            header = str(record.description or "")
            parts = header.split("|") if header else []

            def part(i: int) -> str:
                try:
                    return parts[i]
                except Exception:
                    return ""

            accession = part(0)
            contig = part(1)
            location = part(2)
            strand = part(3)
            protein_id = part(4)
            description = part(5)
            sequence = str(record.seq)

            writer.writerow([
                accession,
                protein_id,
                contig,
                location,
                strand,
                description,
                sequence,
            ])

    print(f"Done! Extracted protein sequences to {output_tsv}")


if __name__ == "__main__":
    main()
