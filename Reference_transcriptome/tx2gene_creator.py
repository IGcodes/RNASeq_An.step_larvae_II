#!/usr/bin/env python3
import csv
import re
import sys

# ---- CONFIG ----
FASTA_FILE = "Anstep_UCI_V1.0_rna.fna"
CSV_OUT    = "Anstep_UCI_V1.0_transcripts_info.csv"

# regex to pull out the gene locus ID in parentheses, e.g. (LOC118517740)
GENE_RE = re.compile(r'\((LOC[0-9A-Za-z_]+)\)')

# ---- PARSING ----
def parse_header(header):
    """
    Given a FASTA header (without the leading '>'),
    extract:
      - accession
      - predicted (bool)
      - organism (first two words)
      - product (text up to ' (LOC...')
      - gene_id (LOC...)
      - molecule (last comma-separated token)
    """
    # Split off accession
    try:
        accession, rest = header.split(maxsplit=1)
    except ValueError:
        # no space after accession?
        return {
            "accession": header,
            "predicted": "",
            "organism": "",
            "product": "",
            "gene_id": "",
            "molecule": ""
        }

    # Detect “PREDICTED:” prefix
    predicted = False
    if rest.startswith("PREDICTED:"):
        predicted = True
        rest = rest[len("PREDICTED:"):].strip()

    # Grab organism as first two words
    parts = rest.split(maxsplit=2)
    organism = " ".join(parts[:2]) if len(parts) >= 2 else parts[0]
    remainder = parts[2] if len(parts) == 3 else ""

    # Pull out gene_id
    m = GENE_RE.search(remainder)
    gene_id = m.group(1) if m else ""

    # Product is everything before " (LOC"
    product = remainder.split(" (LOC", 1)[0].strip()

    # Molecule type: assume it's the last comma-separated field
    molecule = remainder.split(",")[-1].strip()

    return {
        "accession": accession,
        "predicted":  "yes" if predicted else "no",
        "organism":   organism,
        "product":    product,
        "gene_id":    gene_id,
        "molecule":   molecule
    }

# ---- MAIN ----
def main():
    with open(FASTA_FILE, "r") as fasta, \
         open(CSV_OUT, "w", newline="") as csvfile:

        fieldnames = ["accession","predicted","organism","product","gene_id","molecule"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for line in fasta:
            if not line.startswith(">"):
                continue
            header = line[1:].rstrip("\n")
            info = parse_header(header)
            writer.writerow(info)

    print(f"Wrote header info for each transcript to '{CSV_OUT}'")

if __name__ == "__main__":
    main()
