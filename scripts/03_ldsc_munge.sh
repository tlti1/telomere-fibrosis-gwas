#!/bin/bash
# =============================================================================
# Script: 03_ldsc_munge.sh
# Author: Teresa Ironman
# Date: 2026
# Description: Reformats GWAS summary statistics into LDSC format using
#              munge_sumstats.py. Must be run within the ldsc conda environment
#              (Python 2.7). Processes telomere and all 12 fibrosis datasets.
#              Sample sizes per organ from Joof et al. 2026 Table 1.
# Usage: conda activate ldsc && bash 03_ldsc_munge.sh
# =============================================================================

LDSC_DIR="/scratch/alice/t/tlti1/ldsc"
TELOMERE_QC="/scratch/alice/t/tlti1/gwas/telomere/telomere_QC.txt"
FIBROSIS_DIR="/scratch/alice/t/tlti1/gwas/12organ"
MUNGED_DIR="/scratch/alice/t/tlti1/gwas/munged"
HM3="/scratch/alice/t/tlti1/hm3_no_MHC.list.txt"
mkdir -p "$MUNGED_DIR"

declare -A N
N[bile]=78586; N[cardiovascular]=370381; N[diabetes]=356334
N[integumentary]=125148; N[intestinalpanc]=419221; N[liver]=79886
N[lymphatic]=3674; N[reproductive]=85943; N[respiratory]=53167
N[skeletal]=168849; N[systemic]=46630; N[urinary]=546819

echo "Munging telomere summary statistics..."
python "$LDSC_DIR/munge_sumstats.py" \
    --sumstats "$TELOMERE_QC" \
    --out "/scratch/alice/t/tlti1/gwas/telomere/telomere_munged" \
    --snp rs_id --a1 effect_allele --a2 other_allele \
    --p p_value --frq effect_allele_frequency --n-col n

echo "Munging fibrosis summary statistics..."
for organ in bile cardiovascular diabetes integumentary intestinalpanc \
             liver lymphatic reproductive respiratory skeletal systemic urinary; do
    echo "  $organ (N=${N[$organ]})"
    python "$LDSC_DIR/munge_sumstats.py" \
        --sumstats "$FIBROSIS_DIR/results_${organ}_european_info0.8.txt" \
        --out "$MUNGED_DIR/${organ}_munged" \
        --snp rsid --a1 effect_allele --a2 reference_allele \
        --p p --frq eaf --N "${N[$organ]}"
done
echo "Munging complete."
