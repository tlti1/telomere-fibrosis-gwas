#!/bin/bash
# =============================================================================
# Script: 04_ldsc_rg.sh
# Author: Teresa Ironman
# Date: 2026
# Description: Estimates genome-wide genetic correlation between telomere
#              length and fibrosis across 12 organ systems using LDSC v1.0.1.
#              The --no-check-alleles flag is applied due to insertion-deletion
#              variants causing allele incompatibilities between munged files.
#              Bonferroni correction threshold: p < 0.05/12 = 0.0042.
# Usage: conda activate ldsc && bash 04_ldsc_rg.sh
# =============================================================================

LDSC_DIR="/scratch/alice/t/tlti1/ldsc"
TELOMERE="/scratch/alice/t/tlti1/gwas/telomere/telomere_munged.sumstats.gz"
MUNGED_DIR="/scratch/alice/t/tlti1/gwas/munged"
LDSCORE="/scratch/alice/t/tlti1/LDscore/LDscore."
RESULTS_DIR="/scratch/alice/t/tlti1/gwas/ldsc_results"
mkdir -p "$RESULTS_DIR"

echo "Running LDSC genetic correlation analyses..."
for organ in bile cardiovascular diabetes integumentary intestinalpanc \
             liver lymphatic reproductive respiratory skeletal systemic urinary; do
    echo "Processing: $organ"
    python "$LDSC_DIR/ldsc.py" \
        --rg "$TELOMERE,$MUNGED_DIR/${organ}_munged.sumstats.gz" \
        --ref-ld-chr "$LDSCORE" --w-ld-chr "$LDSCORE" \
        --no-check-alleles \
        --out "$RESULTS_DIR/${organ}_rg"
done
echo "All analyses complete."
