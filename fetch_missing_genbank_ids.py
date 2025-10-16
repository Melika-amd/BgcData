"""
fetch_missing_genbank_ids.py (robust version)

✔ Handles any input ID (GenBank, RefSeq, UniProt, locus tags, or gene names)
✔ Tries both [Accession] and [All Fields]
✔ Retries automatically when NCBI closes the connection
✔ Writes the resolved GenBank accession version into a new TSV column
"""

import csv
import time
import random
from pathlib import Path
from typing import Dict, Optional

import requests
from requests import Response, Session



SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_FILE = SCRIPT_DIR / "protein_id_types_with_accession.tsv"
OUTPUT_FILE = SCRIPT_DIR / "unique_proteins_with_genbank.tsv"
EMAIL = "melika@example.com"          
API_KEY = None                        
DELAY = 0.34                          
MAX_RETRIES = 3                       



EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def entrez_get(endpoint: str, params: dict, session: Session, timeout: int = 30) -> Response:
    """Perform GET request to NCBI E-utilities endpoint with automatic retry and backoff."""
    q = dict(params)
    q["email"] = EMAIL
    if API_KEY:
        q["api_key"] = API_KEY
    url = f"{EUTILS_BASE}/{endpoint}"

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            r = session.get(url, params=q, timeout=timeout)
            r.raise_for_status()
            return r
        except (requests.exceptions.ConnectionError, requests.exceptions.Timeout) as e:
            print(f" Connection error (attempt {attempt}/{MAX_RETRIES}) — {e}")
            if attempt < MAX_RETRIES:
                sleep_time = DELAY * 2 + random.uniform(0.5, 2.0)
                print(f"   🔄 Retrying in {sleep_time:.1f}s...")
                time.sleep(sleep_time)
            else:
                print(f" Giving up on request after {MAX_RETRIES} failed attempts.")
                raise
        except requests.exceptions.HTTPError as e:
            print(f" HTTP error for {url}: {e}")
            raise


def esearch_protein_id(protein_id: str, session: Session) -> Optional[str]:
    """
    Return Entrez UID for a given protein ID.
    Tries [Accession] first, then falls back to [All Fields].
    """
    import xml.etree.ElementTree as ET

    for field in ["Accession", "All Fields"]:
        params = {
            "db": "protein",
            "retmax": "1",
            "term": f"{protein_id}[{field}]",
        }
        try:
            r = entrez_get("esearch.fcgi", params, session)
            root = ET.fromstring(r.text)
            idlist = root.find("IdList")
            if idlist is not None:
                for elem in idlist.findall("Id"):
                    if elem.text:
                        return elem.text.strip()
        except Exception as e:
            print(f" esearch failed for {protein_id} ({field}): {e}")
            continue
        if DELAY:
            time.sleep(DELAY)
    return None


def esummary_accession_version(uid: str, session: Session) -> Optional[str]:
    """Return AccessionVersion for a given UID."""
    params = {"db": "protein", "id": uid, "retmode": "json"}
    try:
        r = entrez_get("esummary.fcgi", params, session)
        data = r.json()
        record = data.get("result", {}).get(uid)
        if not record:
            return None
        return record.get("accessionversion")
    except Exception as e:
        print(f"  esummary failed for UID {uid}: {e}")
        return None


def resolve_genbank_id(pid: str, session: Session, cache: Dict[str, Optional[str]]) -> Optional[str]:
    """Find GenBank accession for one protein/locus ID, using cache to avoid duplicates."""
    if pid in cache:
        return cache[pid]

    uid = esearch_protein_id(pid, session)
    if not uid:
        cache[pid] = None
        return None

    if DELAY:
        time.sleep(DELAY)

    acc = esummary_accession_version(uid, session)
    if DELAY:
        time.sleep(DELAY)

    cache[pid] = acc
    return acc


def main():
    if not INPUT_FILE.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_FILE}")

    print(f"🔍 Reading {INPUT_FILE}")
    session = requests.Session()
    cache: Dict[str, Optional[str]] = {}

    with INPUT_FILE.open("r", newline="", encoding="utf-8") as infile, OUTPUT_FILE.open(
        "w", newline="", encoding="utf-8"
    ) as outfile:
        reader = csv.DictReader(infile, delimiter="\t")
        if not {"accession", "protein_id", "datatype"}.issubset(reader.fieldnames or []):
            raise ValueError("Input must contain 'accession', 'protein_id', and 'datatype' columns.")

        fieldnames = list(reader.fieldnames) + ["new_genbank_id"]
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        total = 0
        for row in reader:
            total += 1
            pid = row.get("protein_id", "")
            datatype = row.get("datatype", "")
            new_id = ""

            if pid and "not in ncbi" in datatype.lower():
                try:
                    new_id = resolve_genbank_id(pid, session, cache) or ""
                    print(f"{total:05d} | {pid:<20} → {new_id or 'Not found'}")
                except Exception as e:
                    print(f" Skipping {pid} due to error: {e}")
                    new_id = ""
            else:
                new_id = ""

            row["new_genbank_id"] = new_id
            writer.writerow(row)

    print(f"\nDone! Results saved to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
