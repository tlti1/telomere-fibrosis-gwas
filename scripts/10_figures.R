#!/usr/bin/env Rscript
# ============================================================================
# Script: 10_figures.R
# Author: Teresa Ironman
# Date: 2026
# Description: Produces two summary figures:
#   1. Bar chart of direction of effect consistency per organ system
#   2. Forest plot of LDSC genetic correlations across 12 organ systems
# ============================================================================

library(ggplot2)

OUT_DIR <- "/home/tlti1/DISS/figures"
dir.create(OUT_DIR, showWarnings=FALSE)

# ---- Figure 1: Direction of effect bar chart --------------------------------
dir_data <- data.frame(
  organ = c("Bile","Cardiovascular","Diabetes","Integumentary",
            "Intestinal-Pancreatic","Liver","Lymphatic","Reproductive",
            "Respiratory","Skeletal","Systemic","Urinary"),
  pct = c(49.6, 59.6, 43.8, 51.0, 52.7, 51.0,
          52.2, 50.2, 63.5, 58.6, 53.7, 55.5)
)

p1 <- ggplot(dir_data, aes(x=reorder(organ, pct), y=pct)) +
  geom_col(fill="#2c7bb6", alpha=0.8) +
  geom_hline(yintercept=50, linetype="dashed", colour="red", linewidth=0.6) +
  coord_flip() +
  labs(title="Direction of Effect Consistency Between Telomere Length and Fibrosis",
       subtitle="Dashed line = 50% (chance level)",
       x="Organ System", y="Variants with Consistent Direction (%)") +
  theme_minimal() +
  theme(plot.title=element_text(face="bold", size=11))

png(file.path(OUT_DIR, "direction_of_effect_barchart.png"),
    width=800, height=500, res=120)
print(p1)
dev.off()
cat("Saved: direction_of_effect_barchart.png
")

# ---- Figure 2: LDSC forest plot --------------------------------------------
ldsc <- data.frame(
  organ = c("Bile","Cardiovascular","Diabetes","Integumentary",
            "Intestinal-Pancreatic","Liver","Reproductive",
            "Respiratory","Skeletal","Systemic","Urinary"),
  rg = c(0.0107,-0.0604,-0.0266,-0.0016,-0.0219,0.0013,
         0.0128,-0.1868,-0.0719,-0.0548,-0.0728),
  se = c(0.0355,0.0181,0.0155,0.0314,0.0164,0.0383,
         0.0351,0.072,0.0324,0.0467,0.0392),
  p  = c(0.763,0.0008,0.087,0.959,0.181,0.974,
         0.714,0.0095,0.026,0.241,0.063)
)
ldsc$ci_lower <- ldsc$rg - 1.96 * ldsc$se
ldsc$ci_upper <- ldsc$rg + 1.96 * ldsc$se
ldsc$sig <- ifelse(ldsc$p < 0.0042, "Bonferroni significant",
            ifelse(ldsc$p < 0.05,   "Nominally significant", "Not significant"))
ldsc$sig <- factor(ldsc$sig,
  levels=c("Bonferroni significant","Nominally significant","Not significant"))

p2 <- ggplot(ldsc, aes(x=rg, y=reorder(organ, rg), colour=sig)) +
  geom_vline(xintercept=0, linetype="dashed", colour="grey50") +
  geom_errorbarh(aes(xmin=ci_lower, xmax=ci_upper), height=0.3) +
  geom_point(size=3) +
  scale_colour_manual(values=c("Bonferroni significant"="#c0392b",
                                "Nominally significant"="#e67e22",
                                "Not significant"="grey50")) +
  labs(title="Genetic Correlation Between Telomere Length and Fibrosis",
       subtitle="Estimated using LD Score Regression (LDSC)",
       x="Genetic Correlation (rg)", y="Organ System", colour="Significance") +
  theme_minimal() +
  theme(plot.title=element_text(face="bold", size=12),
        legend.position="bottom")

png(file.path(OUT_DIR, "ldsc_forest_plot.png"),
    width=800, height=500, res=120)
print(p2)
dev.off()
cat("Saved: ldsc_forest_plot.png
")
cat("All figures complete!
")
