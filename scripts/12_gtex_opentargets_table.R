#!/usr/bin/env Rscript
# =============================================================================
# Script: 12_gtex_opentargets_table.R
# Author: Teresa Ironman
# Date: 2026
# Description: Compiles bioinformatic follow-up results from GTEx (V10) and
#              Open Targets Genetics for the 6 strongly colocalised regions
#              (PP.H4 >= 0.8). Includes eQTL/sQTL evidence, likely causal
#              gene assignments based on allele concordance checks, V2G scores
#              from Open Targets, and biological mechanism summaries.
#              Allele concordance was verified by confirming that the GTEx
#              normalised effect size (NES) direction was consistent with the
#              telomere and fibrosis effect directions in our datasets.
#              Input: Manual GTEx and Open Targets lookups
#              Output: gtex_followup_table_updated.xlsx, opentargets_table.xlsx
# =============================================================================

library(writexl)

OUTPUT_DIR <- "/home/tlti1/DISS/coloc"

# GTEx follow-up table
gtex_table <- data.frame(
  Organ = c("Respiratory","Cardiovascular","Cardiovascular",
            "Diabetes","Diabetes","Diabetes"),
  Lead_variant = c("rs7705526","rs11556924","rs4766578",
                   "rs4766578","rs3768321","rs35601737"),
  Annotated_gene = c("TERT","ZC3HC1","ATXN2","ATXN2","PABPC4","TRMT1"),
  PP.H4 = c(0.816, 0.987, 0.970, 0.893, 0.992, 0.917),
  Likely_causal_gene = c("TERT","ZC3HC1","ATXN2","ATXN2","PABPC4","TRMT1"),
  GTEx_evidence = c(
    "No TERT eQTL in GTEx; literature confirms A allele increases TERT enhancer activity",
    "T allele decreases ZC3HC1 in fibroblasts (NES=-0.066); KLHDC10 stronger in lung but less tissue-relevant",
    "T allele decreases ATXN2 (NES=+0.066 for A allele); ALDH2 excluded - inconsistent allele direction",
    "T allele decreases ATXN2 (NES=+0.066 for A allele); ALDH2 excluded - inconsistent allele direction",
    "T allele decreases PABPC4 in fibroblasts (NES=-0.43, p=7e-44) and multiple tissues",
    "G allele increases TRMT1 splicing across all tissues (sQTL NES=+1.3, p=1e-317)"
  ),
  Mechanism = c(
    "A allele increases TERT enhancer activity -> longer telomeres -> less fibrosis",
    "T allele decreases ZC3HC1 in fibroblasts -> longer telomeres -> less cardiovascular fibrosis",
    "T allele decreases ATXN2 -> shorter telomeres -> more cardiovascular fibrosis",
    "T allele decreases ATXN2 -> shorter telomeres -> more diabetes fibrosis; same locus as SH2B3 (Joof et al. 2026)",
    "T allele decreases PABPC4 in fibroblasts and artery -> shorter telomeres -> more diabetes fibrosis",
    "G allele alters TRMT1 splicing -> impaired tRNA modification -> shorter telomeres -> more diabetes fibrosis"
  ),
  Allele_concordant = c("Yes","Yes",
                        "Yes - after allele check (not ALDH2)",
                        "Yes - after allele check (not ALDH2)",
                        "Yes","Likely yes")
)

write_xlsx(gtex_table, file.path(OUTPUT_DIR, "gtex_followup_table_updated.xlsx"))
cat("Saved: gtex_followup_table_updated.xlsx\n")

# Open Targets table
ot_table <- data.frame(
  Organ = c("Respiratory","Cardiovascular","Cardiovascular & Diabetes",
            "Diabetes","Diabetes"),
  Lead_variant = c("rs7705526","rs11556924","rs4766578",
                   "rs3768321","rs35601737"),
  PP.H4 = c(0.816, 0.987, 0.970, 0.992, 0.917),
  Top_L2G_gene = c("TERT","ZC3HC1","SH2B3","PABPC4","TRMT1"),
  L2G_score = c(0.818, 0.931, 0.984, 0.893, 0.677),
  Variant_consequence = c(
    "Intron variant in TERT - enhancer activity confirmed",
    "Missense variant in ZC3HC1 - R363H amino acid change",
    "Intron variant in ATXN2 - regulatory",
    "Intron variant near PABPC4",
    "Intron variant in TRMT1 - specific transcript isoform affected"
  ),
  Key_GWAS_associations = c(
    "Dyskeratosis congenita; Idiopathic pulmonary fibrosis (ClinVar)",
    "Ischemic heart disease; Coronary artery bypass; Blood pressure",
    "Cardiovascular disorder; Hypertension; Rheumatoid arthritis",
    "Type 2 diabetes; C-reactive protein; HDL cholesterol",
    "Urea levels; eGFR; Platelet count"
  ),
  Key_molQTL = c(
    "Enhancer predictions for TERT across multiple tissues",
    "eQTL for ZC3HC1 in thyroid (credible set size=1)",
    "SH2B3 L2G=0.984; pQTL for immune proteins in blood",
    "eQTL for PABPC4 in fibroblasts (credible set size=1); eQTL in pancreas",
    "Transcript QTL for TRMT1 across all tissues (p=1e-181, credible set=1)"
  ),
  Notes = c(
    "ClinVar confirms IPF and dyskeratosis congenita associations",
    "Missense R363H directly changes ZC3HC1 protein - strongest functional evidence",
    "Open Targets points to SH2B3 not ATXN2/ALDH2; converges with Joof et al. 2026",
    "Cleanest result - credible set size 1 in fibroblasts and pancreas",
    "Lower L2G score (0.677); enhancer predictions also suggest LYL1 as alternative"
  )
)

write_xlsx(ot_table, file.path(OUTPUT_DIR, "opentargets_table.xlsx"))
cat("Saved: opentargets_table.xlsx\n")
cat("Bioinformatic follow-up tables complete!\n")
