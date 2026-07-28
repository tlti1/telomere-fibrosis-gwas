#!/bin/bash
# =============================================================================
# Script: 06_coloc_regions_extract.sh
# Author: Teresa Ironman
# Date: 2026
# Description: Extracts +-500kb windows around each Bonferroni-significant
#              lead variant from telomere and fibrosis summary statistics
#              for colocalisation analysis. MAF >= 1% filter applied.
#              Input: Telomere QC file, fibrosis summary statistics
#              Output: Region files for coloc analysis
# =============================================================================

TELOMERE_QC="/scratch/alice/t/tlti1/gwas/telomere/telomere_QC.txt"
FIBROSIS_DIR="/scratch/alice/t/tlti1/gwas/12organ"
OUTPUT_DIR="/scratch/alice/t/tlti1/gwas/coloc_regions"
mkdir -p "$OUTPUT_DIR"

# Format: organ:gene:chr:pos:start:end (+-500kb windows)
REGIONS=(
    "respiratory:TERT:5:1285974:785974:1785974"
    "cardiovascular:ZC3HC1:7:129663496:129163496:130163496"
    "cardiovascular:ATXN2:12:111904371:111404371:112404371"
    "diabetes:HHEX:10:94446635:93946635:94946635"
    "diabetes:PABPC4:1:40035928:39535928:40535928"
    "diabetes:MST1R:3:49936102:49436102:50436102"
    "diabetes:SEC61A2:10:12201002:11701002:12701002"
    "diabetes:TRMT1:19:13220703:12720703:13720703"
    "diabetes:ATXN2:12:111904371:111404371:112404371"
    "intestinalpanc:DMC1:22:38921084:38421084:39421084"
)

for region in "${REGIONS[@]}"; do
    organ=$(echo "$region" | cut -d: -f1)
    gene=$(echo "$region" | cut -d: -f2)
    chr=$(echo "$region" | cut -d: -f3)
    start=$(echo "$region" | cut -d: -f5)
    end=$(echo "$region" | cut -d: -f6)
    echo "Extracting: $organ $gene (chr${chr}:${start}-${end})"

    # Fibrosis region with MAF >= 1%
    awk -F'\t' -v c="$chr" -v s="$start" -v e="$end" \
        'NR==1||($1==c&&$2>=s&&$2<=e&&$6>=0.01&&$6<=0.99)' \
        "$FIBROSIS_DIR/results_${organ}_european_info0.8.txt" \
        > "$OUTPUT_DIR/${organ}_${gene}.txt"

    # Telomere region with MAF >= 1%
    awk -F'\t' -v c="$chr" -v s="$start" -v e="$end" \
        'NR==1||($1==c&&$2>=s&&$2<=e&&$7>=0.01&&$7<=0.99)' \
        "$TELOMERE_QC" > "$OUTPUT_DIR/telomere_${gene}.txt"
done
echo "Region extraction complete. Files saved to $OUTPUT_DIR"
