#!/bin/bash
# =============================================================================
# Script: 07_vep_annotation.sh
# Author: Teresa Ironman
# Date: 2026
# Description: Extracts rsIDs from Bonferroni-significant variants and
#              submits to Ensembl Variant Effect Predictor (VEP) for gene
#              annotation. VEP was run via the online tool at
#              https://www.ensembl.org/Tools/VEP using Ensembl/GENCODE
#              transcripts. MANE Select transcripts were preferred where
#              available. This script extracts the rsID list and parses
#              the VEP output to assign one gene per variant.
#              Input: Bonferroni-significant variants file
#              Output: rsID list for VEP, parsed gene annotations
# =============================================================================

INPUT="/scratch/alice/t/tlti1/gwas/significant/significant_all_organs.txt"
OUTPUT_DIR="/scratch/alice/t/tlti1/gwas/significant"

echo "Extracting rsIDs for VEP annotation..."

# Extract unique rsIDs from Bonferroni-significant variants
awk -F'\t' 'NR>1 {print $4}' "$INPUT" | sort | uniq \
    > "$OUTPUT_DIR/rsids_for_vep.txt"

echo "$(wc -l < $OUTPUT_DIR/rsids_for_vep.txt) unique rsIDs extracted"
echo "Submit $OUTPUT_DIR/rsids_for_vep.txt to Ensembl VEP online:"
echo "https://www.ensembl.org/Tools/VEP"
echo ""
echo "Settings used:"
echo "  - Transcript database: Ensembl/GENCODE transcripts"
echo "  - Additional annotations: MANE Select"
echo ""
echo "After downloading VEP output, parse with:"
echo "  bash 07b_parse_vep.sh <vep_output_file>"
