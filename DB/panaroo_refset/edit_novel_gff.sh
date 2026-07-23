#!/bin/bash

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 input.gff" >&2
    exit 2
fi

gff="$1"

if [[ ! -f "$gff" ]]; then
    echo "ERROR: file not found: $gff" >&2
    exit 1
fi

if [[ "$gff" != *.gff ]]; then
    echo "ERROR: Panaroo 1.5.2 requires the input filename to end in .gff: $gff" >&2
    exit 1
fi

# Panaroo-generate-gffs treats Bakta's single-# metadata lines as GFF
# feature rows. Convert only the standard Bakta metadata comments to ##.
sed -i -E \
    '/^# (Annotated with Bakta|Software:|Database:|DOI:|URL:)/s/^#/##/' \
    "$gff"

# Add attributes required for Bakta GFF input to Panaroo 1.5.2.
# The guards make this operation idempotent.
sed -i '
    /\tCDS\t/ {
        /;panaroo_id=/! s/;product=/;panaroo_id=;product=/
        /;prepanaroo_inference=/! s/;panaroo_id=;/;panaroo_id=;prepanaroo_inference=Unknown_inference;/
    }
' "$gff"

cds=$(awk -F '\t' '$3 == "CDS" {n++} END {print n+0}' "$gff")
panaroo=$(awk -F '\t' \
    '$3 == "CDS" && $9 ~ /(^|;)panaroo_id=/ {n++} END {print n+0}' \
    "$gff")
inference=$(awk -F '\t' \
    '$3 == "CDS" && $9 ~ /(^|;)prepanaroo_inference=/ {n++} END {print n+0}' \
    "$gff")

if [[ "$cds" -ne "$panaroo" || "$cds" -ne "$inference" ]]; then
    echo "ERROR: incomplete Panaroo attributes in $gff" >&2
    echo "CDS=$cds panaroo_id=$panaroo inference=$inference" >&2
    exit 1
fi

echo "$gff CDS=$cds panaroo_id=$panaroo inference=$inference"
