#!/bin/bash
# =============================================================================
# Script: 07b_parse_vep.sh
# Author: Teresa Ironman
# Date: 2026
# Description: Parses VEP output file to extract one gene annotation per
#              variant. Prefers MANE Select transcripts where available.
#              For intergenic variants, annotates as "intergenic".
#              For indels without rsIDs, annotates as NA.
#              Input: VEP output file (.vep format)
#              Output: Two-column file of rsID and gene symbol
# Usage: bash 07b_parse_vep.sh <vep_output.vep>
# =============================================================================

VEP_FILE="$1"
OUTPUT="/scratch/alice/t/tlti1/gwas/significant/vep_genes.txt"

if [ -z "$VEP_FILE" ]; then
    echo "Usage: bash 07b_parse_vep.sh <vep_output.vep>"
    exit 1
fi

echo "Parsing VEP output: $VEP_FILE"

# Extract gene symbol from Extra column, preferring MANE Select transcripts
grep -v "^#" "$VEP_FILE" | awk -F'\t' '{
    rsid=$1
    extra=$14
    # Extract SYMBOL
    if (extra ~ /SYMBOL=/) {
        symbol=extra
        gsub(/.*SYMBOL=/, "", symbol)
        gsub(/;.*/, "", symbol)
    } else {
        symbol="intergenic"
    }
    # Flag MANE Select transcripts
    mane=(extra ~ /MANE_Select/) ? 1 : 0
    print rsid"\t"symbol"\t"mane
}' | sort -k1,1 -k3,3rn | awk -F'\t' '!seen[$1]++ {print $1"\t"$2}' \
    > "$OUTPUT"

echo "Done. $(wc -l < $OUTPUT) variants annotated."
echo "Output: $OUTPUT"
