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
        default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "output", "mibig_fasta_proteins.tsv"),
        help="Path to output TSV file",
    )
    args = parser.parse_args()

    input_fasta = args.input
    output_tsv = args.output

    # Ensure parent directory exists for output
    os.makedirs(os.path.dirname(os.path.abspath(output_tsv)), exist_ok=True)

    with open(output_tsv, "w", newline="", encoding="utf-8") as tsv_file:
        writer = csv.writer(tsv_file, delimiter="\t")
        # Fixed header: accession
        writer.writerow(["accession", "protein_id", "description", "sequence"])  

        for record in SeqIO.parse(input_fasta, "fasta"):
            header = record.description
            accession = record.id.split("_prot_")[0]
            protein_id = record.id
            description = header
            sequence = str(record.seq)

            writer.writerow([accession, protein_id, description, sequence])

    print(f"Done! Extracted protein sequences to {output_tsv}")


if __name__ == "__main__":
    main()
