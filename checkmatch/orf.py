import requests
import pandas as pd
import time

def fetch_uniprot_entry(protein_id):
    url = f"https://rest.uniprot.org/uniprotkb/search?query={protein_id}&format=json"
    try:
        r = requests.get(url, timeout=10)
        r.raise_for_status()
        data = r.json()
        if "results" in data and len(data["results"]) > 0:
            return data["results"][0]
    except Exception as e:
        print(f"[ERROR] {protein_id}: {e}")
    return None

def extract_orf_name(entry):
    try:
        for g in entry.get("genes", []):
            if "orfNames" in g and len(g["orfNames"]) > 0:
                return g["orfNames"][0]["value"]
    except Exception:
        pass
    return None

def main():
    
    df = pd.read_csv("checkmatch/final_output/unique_proteins.tsv", sep="\t")

    orf_names = []
    for idx, pid in enumerate(df["protein_id"], start=1):
        entry = fetch_uniprot_entry(pid)
        if entry:
            orf_name = extract_orf_name(entry)
            orf_names.append(orf_name)
            print(f"[{idx}/{len(df)}] {pid} -> ORF: {orf_name}")
        else:
            orf_names.append(None)

        
        time.sleep(0.5)

    df["orf_name"] = orf_names
    df.to_csv("checkmatch/final_output/unique_proteins_with_orf.tsv", sep="\t", index=False)
    print("[DONE] Results saved to checkmatch/final_output/unique_proteins_with_orf.tsv")

if __name__ == "__main__":
    main()
