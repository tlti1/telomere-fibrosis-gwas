#!/bin/bash
# =============================================================================
# Script: 01_qc_filter.sh
# Author: Teresa Ironman
# Date: 2026
# Description: Quality control filtering of telomere length GWAS summary 
#              statistics. Filters variants by MAF >= 1% and INFO > 0.5.
#              Input: Raw telomere GWAS summary statistics (GCST90435144)
#              Output: QC-filtered summary statistics
# =============================================================================

TELOMERE_RAW="/scratch/alice/t/tlti1/gwas/telomere/GCST90435144.tsv.gz"
TELOMERE_QC="/scratch/alice/t/tlti1/gwas/telomere/telomere_QC.txt"

echo "Starting QC filtering of telomere summary statistics..."
zcat "$TELOMERE_RAW" | awk -F'\t' '
NR==1 {print; next}
$7 >= 0.01 && $7 <= 0.99 && $9 > 0.5 {print}
' > "$TELOMERE_QC"
echo "Done. Variants remaining: $(wc -l < $TELOMERE_QC)"
