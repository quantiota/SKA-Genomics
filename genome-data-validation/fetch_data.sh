#!/usr/bin/env bash
# Download the E. coli K-12 MG1655 reference genome (NC_000913.3) from NCBI.
set -euo pipefail
mkdir -p data
OUT="data/GCF_000005845.2_ASM584v2_genomic.fna"
echo "Fetching NC_000913.3 from NCBI ..."
curl -sS -o "$OUT" \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=fasta&retmode=text"
echo "Wrote $OUT  ($(grep -vc '^>' "$OUT") sequence lines)"
head -1 "$OUT"
