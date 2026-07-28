#!/usr/bin/env Rscript
# =============================================================================
# Script: 11_bonferroni_table.R
# Author: Teresa Ironman
# Date: 2026
# Description: Builds the final annotated table of Bonferroni-significant
#              overlapping variants between telomere length and fibrosis.
#              Adds VEP gene annotations, confidence intervals for both
#              fibrosis and telomere effect sizes, and allele status.
#              Both fibrosis and telomere betas reported for the same
#              effect allele (telomere beta sign corrected where alleles
#              were flipped between datasets).
#              Input: Overlap corrected files, VEP gene annotations
#              Output: bonferroni_table_final.xlsx
# =============================================================================

library(writexl)

OVERLAP_DIR <- "/scratch/alice/t/tlti1/gwas/overlap_corrected"
VEP_FILE    <- "/scratch/alice/t/tlti1/gwas/significant/vep_genes.txt"
OUTPUT      <- "/home/tlti1/DISS/coloc/bonferroni_table_final.xlsx"
THRESHOLD   <- 0.0001262626  # Bonferroni: 0.05 / 396

# Read VEP annotations
vep <- read.table(VEP_FILE, header=FALSE, sep="\t",
                  col.names=c("rsid", "gene"))

# Read all organ overlap files and filter to Bonferroni threshold
all_results <- data.frame()

organs <- c("bile","cardiovascular","diabetes","integumentary","intestinalpanc",
            "liver","lymphatic","reproductive","respiratory","skeletal",
            "systemic","urinary")

for (organ in organs) {
  f <- file.path(OVERLAP_DIR,
                 paste0("results_", organ, "_european_info0.8_corrected.txt"))
  if (!file.exists(f)) next
  
  dat <- read.table(f, header=TRUE, sep="\t", fill=TRUE)
  dat <- dat[!is.na(dat$fibrosis_p) & dat$fibrosis_p < THRESHOLD, ]
  
  if (nrow(dat) == 0) next
  dat$organ <- organ
  all_results <- rbind(all_results, dat)
}

# Add confidence intervals (beta +/- 1.96 * SE)
all_results$fibrosis_ci_lower  <- round(all_results$fibrosis_beta  - 
                                         1.96 * all_results$fibrosis_se, 4)
all_results$fibrosis_ci_upper  <- round(all_results$fibrosis_beta  + 
                                         1.96 * all_results$fibrosis_se, 4)
all_results$telomere_ci_lower  <- round(all_results$telomere_beta  - 
                                         1.96 * all_results$telomere_se, 4)
all_results$telomere_ci_upper  <- round(all_results$telomere_beta  + 
                                         1.96 * all_results$telomere_se, 4)
all_results$fibrosis_beta      <- round(all_results$fibrosis_beta, 4)
all_results$telomere_beta      <- round(all_results$telomere_beta, 4)

# Merge VEP gene annotations
all_results <- merge(all_results, vep, by="rsid", all.x=TRUE)

# Select and order final columns
final <- all_results[, c("organ","chr","position","rsid","gene",
                          "effect_allele",
                          "fibrosis_beta","fibrosis_ci_lower",
                          "fibrosis_ci_upper","fibrosis_p",
                          "telomere_beta","telomere_ci_lower",
                          "telomere_ci_upper","telomere_p",
                          "allele_status","direction")]

final <- final[order(final$organ, final$fibrosis_p), ]

cat("Bonferroni-significant variants:", nrow(final), "\n")
print(final)

write_xlsx(final, OUTPUT)
cat("Saved to:", OUTPUT, "\n")
