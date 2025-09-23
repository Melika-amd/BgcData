Simple Overview
===============

Parsers
-------
- FASTA → TSV (`MIBiG_BGC/FASTA/Fasta_parser.py`)
  - Reads: `MIBiG_BGC/FASTA/FASTAFiles/mibig_prot_seqs_4.0.fasta`
  - Writes: `MIBiG_BGC/output/mibig_fasta_proteins.tsv`
  - Extracts per-protein rows with columns: `accession, protein_id, description, sequence`.

- JSON → TSV (`MIBiG_BGC/JSON/json_to_tsv.py`)
  - Reads: `MIBiG_BGC/JSON/json_files/mibig_json_4.0/`
  - Writes: `MIBiG_BGC/output/mibig_json_metadata.tsv`
  - Aggregates BGC metadata and sorts rows by accession number.

- DeepBGC PFAM cleaner (`DeepBGC/Scripts/DeepBGC_parser.py`)
  - Reads: `DeepBGC/Data/MIBiG.pfam.tsv`
  - Writes: `DeepBGC/Data/cleaned_deepbgc.tsv`
  - Keeps PFAM columns, renames `sequence_id` → `accession`, removes version suffix (e.g., `.1`).

Checkmatch
----------
- Unique finder and merger (`checkmatch/matchfinder.py`)
  - Compares DeepBGC cleaned accessions to FASTA accessions.
  - Outputs to `checkmatch/final_output/`:
    - `unique_deepbgc_accessions.tsv`: BGCs found only in DeepBGC.
    - `combined_unique_plus_mibig.tsv`: union of MIBiG + DeepBGC-unique with a `source` column.

