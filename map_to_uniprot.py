import pandas as pd
import requests
import time


INPUT_FILE = "unique_proteins.tsv"   
OUTPUT_FILE = "protein_id_types.tsv"
EMAIL = "melika@example.com"        
DELAY = 0.34                       



def classify_by_prefix(pid):
    """Quick initial classification based on prefix pattern."""
    if pid.startswith(("WP_", "NP_", "XP_", "YP_")):
        return "RefSeq"
    if len(pid) in [6, 10] and pid[0].isalpha() and pid[1].isdigit():
        return "UniProt"
    if pid[:3].isalpha() and pid.count(".") == 1:
        return "GenBank_or_EMBL"
    return "Unknown"


def check_ncbi(pid):
    """Query NCBI ESummary to confirm if it’s an NCBI record and identify type."""
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {"db": "protein", "id": pid, "retmode": "json", "email": EMAIL}

    try:
        r = requests.get(url, params=params, timeout=15)
        r.raise_for_status()
        data = r.json()
        uids = data.get("result", {}).get("uids", [])
        if not uids:
            return "Not in NCBI"

        uid = uids[0]
        rec = data["result"][uid]
        acc = rec.get("accessionversion", "")
        src = rec.get("source", "")

        if acc.startswith(("WP_", "NP_", "XP_", "YP_")):
            return "RefSeq"
        elif "GenBank" in src or acc[:3].isalpha():
            return "GenBank"
        else:
            return "Other (NCBI)"

    except Exception:
        return "Error"


def main():
    df = pd.read_csv(INPUT_FILE, sep="\t")
    ids = df["protein_id"].dropna().astype(str).tolist()

    results = []
    print(f"🔍 Checking {len(ids)} IDs...")

    for i, pid in enumerate(ids, 1):
        prefix_type = classify_by_prefix(pid)
        if prefix_type in ["RefSeq", "UniProt"]:
            final_type = prefix_type
        else:
            final_type = check_ncbi(pid)

        results.append({"protein_id": pid, "id_type": final_type})
        print(f"{i:05d}/{len(ids)} | {pid} → {final_type}")
        time.sleep(DELAY)

    out = pd.DataFrame(results)
    out.to_csv(OUTPUT_FILE, sep="\t", index=False)
    print(f"\n Done! Results saved to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
