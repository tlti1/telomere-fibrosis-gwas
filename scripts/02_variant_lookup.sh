#!/bin/bash
# =============================================================================
# Script: 02_variant_lookup.sh
# Author: Teresa Ironman
# Date: 2026
# Description: Looks up 396 independent telomere signals (Burren et al. 2024
#              Supplementary Table 6) in fibrosis GWAS across 12 organ systems.
#              Corrects allele alignment and calculates direction of effect.
#              Consistent direction = shorter telomeres, more fibrosis risk.
#              Input: Independent telomere signals, fibrosis summary statistics
#              Output: Corrected overlap files per organ system
# =============================================================================

TELO_SIGNALS="/scratch/alice/t/tlti1/gwas/telomere/telomere_independent_PC1_clean.txt"
FIBROSIS_DIR="/scratch/alice/t/tlti1/gwas/12organ"
OUTPUT_DIR="/scratch/alice/t/tlti1/gwas/overlap_corrected"
mkdir -p "$OUTPUT_DIR"

for organ in bile cardiovascular diabetes integumentary intestinalpanc \
             liver lymphatic reproductive respiratory skeletal systemic urinary; do
    echo "Processing: $organ"
    FIBROSIS_FILE="$FIBROSIS_DIR/results_${organ}_european_info0.8.txt"
    OUTPUT_FILE="$OUTPUT_DIR/results_${organ}_european_info0.8_corrected.txt"

    echo -e "chr\tposition\trsid\treference_allele\teffect_allele\t\
fibrosis_eaf\tfibrosis_info\tfibrosis_beta\tfibrosis_se\tfibrosis_p\t\
fibrosis_hwe_p\ttelomere_beta\ttelomere_se\ttelomere_p\tallele_status\tdirection" \
    > "$OUTPUT_FILE"

    awk -F'\t' '
    NR==1{next}
    FNR==NR{
        rsid=$4; telo_ref[$rsid]=$6; telo_alt[$rsid]=$7
        telo_beta[$rsid]=$8; telo_se[$rsid]=$9; telo_p[$rsid]=$10; next
    }
    FNR==1{next}
    $3 in telo_ref {
        rsid=$3; fib_beta=$8; tb=telo_beta[rsid]
        if (telo_alt[rsid]==$5 && telo_ref[rsid]==$4) {
            status="aligned"; corrected_tb=tb
        } else if (telo_alt[rsid]==$4 && telo_ref[rsid]==$5) {
            status="flipped"; corrected_tb=-tb
        } else {status="unclear"; corrected_tb=tb}
        direction=((fib_beta>0 && corrected_tb<0)||(fib_beta<0 && corrected_tb>0)) \
                  ? "consistent" : "inconsistent"
        print $1"\t"$2"\t"rsid"\t"$4"\t"$5"\t"$6"\t"$7"\t"fib_beta"\t"$9"\t"\
              $10"\t"$11"\t"corrected_tb"\t"telo_se[rsid]"\t"telo_p[rsid]"\t"\
              status"\t"direction
    }' "$TELO_SIGNALS" "$FIBROSIS_FILE" >> "$OUTPUT_FILE"
    echo "  Done: $(wc -l < $OUTPUT_FILE) variants"
done
echo "Variant lookup complete."
