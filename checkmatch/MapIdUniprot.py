import requests
import pandas as pd
import time

INPUT_FILE = "checkmatch/final_output/unique_proteins.tsv"   
OUTPUT_FILE = "checkmatch/unique_proteins_with_embl.tsv"

def fetch_uniprot_and_embl(protein_id):
    url = f"https://rest.uniprot.org/uniprotkb/search?query={protein_id}&format=json&fields=accession,xref_embl"
    try:
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        data = r.json()

        if "results" not in data or len(data["results"]) == 0:
            return None, None

        result = data["results"][0]
        uniprot_id = result.get("primaryAccession", "")

        embl_ids = []
        for xref in result.get("uniProtKBCrossReferences", []):
            if xref.get("database") == "EMBL":
                embl_ids.append(xref.get("id"))

        embl_id = ",".join(embl_ids) if embl_ids else None
        return uniprot_id, embl_id

    except Exception as e:
        print(f"[ERROR] {protein_id}: {e}")
        return None, None

def main():
    df = pd.read_csv(INPUT_FILE, sep="\t", dtype=str)
    if "protein_id" not in df.columns:
        raise ValueError("Input file must contain a 'protein_id' column")

    uniprot_ids, embl_ids = [], []

    for i, pid in enumerate(df["protein_id"], 1):
        print(f"[LOG] Processing {i}/{len(df)}: {pid}")
        uniprot_id, embl_id = fetch_uniprot_and_embl(pid)
        uniprot_ids.append(uniprot_id)
        embl_ids.append(embl_id)
        time.sleep(0.5) 

    df["uniprot_id"] = uniprot_ids
    df["embl_id"] = embl_ids

    df.to_csv(OUTPUT_FILE, sep="\t", index=False)
    print(f"[DONE] Saved results to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()
