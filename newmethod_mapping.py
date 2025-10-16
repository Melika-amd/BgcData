"""
fetch_uniprot_from_ncbi.py

Given a TSV with a 'protein_id' column, this script queries NCBI's E-utilities API
to find corresponding UniProtKB accessions (cross-references) for each protein.

Output:
    unique_proteins_with_uniprot.tsv — includes new column 'uniprot_id'
"""

import pandas as pd
import requests
import time
import xml.etree.ElementTree as ET
from pathlib import Path


INPUT_FILE = Path("unique_proteins_normalized.tsv")
OUTPUT_FILE = Path("unique_proteins_with_uniprot.tsv")
EMAIL = "melika@example.com"        
DELAY = 0.34                        


def get_uniprot_from_ncbi(protein_id: str) -> str:
    """
    Query NCBI for a protein accession and return the UniProt cross-reference, if any.
    """
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {"db": "protein", "id": protein_id, "retmode": "xml", "email": EMAIL}

    try:
        r = requests.get(base, params=params, timeout=20)
        r.raise_for_status()
    except requests.RequestException:
        return ""

 
    try:
        root = ET.fromstring(r.text)
        for docsum in root.findall(".//DocSum"):
            for item in docsum.findall("Item"):
                if item.attrib.get("Name") == "Extra" and "UniProtKB" in item.text:
                    
                    parts = item.text.split(":")
                    if len(parts) == 2:
                        return parts[-1].strip()
        return ""
    except ET.ParseError:
        return ""


def main():
    if not INPUT_FILE.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_FILE}")

    df = pd.read_csv(INPUT_FILE, sep="\t")
    if "protein_id" not in df.columns:
        raise ValueError("Input file must contain a 'protein_id' column.")

    ids = df["protein_id"].dropna().astype(str).tolist()
    print(f"🔍 Checking {len(ids)} protein IDs from {INPUT_FILE.name}")

    results = []
    for i, pid in enumerate(ids, start=1):
        uniprot_id = get_uniprot_from_ncbi(pid)
        results.append(uniprot_id)
        print(f"{i:05d}/{len(ids)} | {pid} → {uniprot_id or 'No UniProt ID'}")
        time.sleep(DELAY)

    df["uniprot_id"] = results
    df.to_csv(OUTPUT_FILE, sep="\t", index=False)
    print(f"\n Done! Results saved to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
