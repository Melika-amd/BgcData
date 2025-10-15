import pandas as pd

df = pd.read_csv("MIBiG_BGC/output/mibig_fasta_proteins_fixed.tsv", sep="\t")
unique_ids = df["protein_id"].dropna().unique()

pd.DataFrame(unique_ids, columns=["protein_id"]).to_csv("unique_proteins.tsv", sep="\t", index=False)
print("Extracted unique protein IDs to unique_proteins.tsv")
