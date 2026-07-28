#!/bin/bash
# =============================================================================
# Script: 05_bonferroni_filter.sh
# Author: Teresa Ironman
# Date: 2026
# Description: Applies Bonferroni correction threshold (p < 0.05/396 =
#              1.26e-4) to identify variants significantly associated with
#              both telomere length and fibrosis. All 396 variants are
#              genome-wide significant for telomere length by definition.
#              Input: Corrected overlap files from 02_variant_lookup.sh
#              Output: Bonferroni-significant variants, combined output file
# =============================================================================

INPUT_DIR="/scratch/alice/t/tlti1/gwas/overlap_corrected"
OUTPUT_DIR="/scratch/alice/t/tlti1/gwas/significant"
THRESHOLD="0.0001262626"  # 0.05 / 396
mkdir -p "$OUTPUT_DIR"

echo "Applying Bonferroni threshold: p < $THRESHOLD (0.05/396)"

echo -e "organ\tchr\tposition\trsid\teffect_allele\tfibrosis_beta\tfibrosis_se\t\
fibrosis_p\ttelomere_beta\ttelomere_se\ttelomere_p\tallele_status\tdirection" \
> "$OUTPUT_DIR/significant_all_organs.txt"

for file in "$INPUT_DIR"/*.txt; do
    organ=$(basename "$file" | sed 's/results_//;s/_european_info0.8_corrected.txt//')
    count=$(awk -F'\t' -v t="$THRESHOLD" 'NR>1 && $10<t' "$file" | wc -l)
    echo "$organ: $count"
    awk -F'\t' -v organ="$organ" -v t="$THRESHOLD" \
        'NR>1 && $10<t {print organ"\t"$0}' "$file" \
        >> "$OUTPUT_DIR/significant_all_organs.txt"
done
echo "Done. Results: $OUTPUT_DIR/significant_all_organs.txt"
