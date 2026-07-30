#!/usr/bin/env Rscript
# ============================================================================
# Script: 10_figures.R
# Author: Teresa Ironman
# Date: 2026
# Description: Produces all summary figures for the project:
#   1. Bar chart of direction of effect consistency per organ
#   2. Forest plot of LDSC genetic correlations across 12 organs
#   3. Bar chart of PP.H4 colocalisation probabilities
#   4. Dot plot of Open Targets L2G scores for colocalised variants
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
cat("Saved: Figure 1 - direction_of_effect_barchart.png
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
cat("Saved: Figure 2 - ldsc_forest_plot.png
")

# ---- Figure 3: PP.H4 bar chart ---------------------------------------------
coloc_data <- data.frame(
  region = c("TERT
(Respiratory)","ZC3HC1
(Cardiovascular)",
             "ATXN2
(Cardiovascular)","ATXN2
(Diabetes)",
             "PABPC4
(Diabetes)","MST1R
(Diabetes)",
             "SEC61A2
(Diabetes)","TRMT1
(Diabetes)",
             "DMC1
(Intestinal-
Pancreatic)"),
  PP.H4 = c(0.816,0.987,0.970,0.893,0.992,0.765,0.0001,0.917,0.727),
  conclusion = c("Colocalises","Colocalises","Colocalises",
                 "Colocalises","Colocalises","Borderline",
                 "H3 - Different variants","Colocalises","Borderline")
)
coloc_data$conclusion <- factor(coloc_data$conclusion,
  levels=c("Colocalises","Borderline","H3 - Different variants"))
coloc_data$region <- factor(coloc_data$region,
  levels=rev(coloc_data$region))

p3 <- ggplot(coloc_data, aes(x=PP.H4, y=region, fill=conclusion)) +
  geom_col(width=0.7) +
  geom_vline(xintercept=0.8, linetype="dashed",
             colour="black", linewidth=0.7) +
  scale_fill_manual(values=c(
    "Colocalises"="#2c7bb6",
    "Borderline"="#fdae61",
    "H3 - Different variants"="#d7191c"
  )) +
  scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.2)) +
  labs(title="Posterior Probability of Shared Causal Variant (PP.H4)",
       subtitle="Dashed line indicates PP.H4 = 0.8 colocalisation threshold",
       x="PP.H4", y="Genomic Region", fill="Conclusion") +
  theme_minimal() +
  theme(plot.title=element_text(face="bold", size=11),
        plot.subtitle=element_text(size=9, colour="grey40"),
        legend.position="bottom",
        axis.text.y=element_text(size=9))

png(file.path(OUT_DIR, "coloc_PPH4_barchart.png"),
    width=800, height=500, res=120)
print(p3)
dev.off()
cat("Saved: Figure 3 - coloc_PPH4_barchart.png
")

# ---- Figure 4: L2G dot plot ------------------------------------------------
l2g_data <- data.frame(
  variant = c("rs7705526","rs11556924","rs4766578","rs3768321","rs35601737"),
  organ   = c("Respiratory","Cardiovascular","Cardiovascular & Diabetes",
              "Diabetes","Diabetes"),
  gene    = c("TERT","ZC3HC1","SH2B3","PABPC4","TRMT1"),
  L2G     = c(0.818, 0.931, 0.984, 0.893, 0.677),
  PP.H4   = c(0.816, 0.987, 0.970, 0.992, 0.917)
)
l2g_data$label <- paste0(l2g_data$gene, "
(", l2g_data$variant, ")")
l2g_data$label <- factor(l2g_data$label,
                          levels=l2g_data$label[order(l2g_data$L2G)])

p4 <- ggplot(l2g_data, aes(x=L2G, y=label, colour=organ, size=PP.H4)) +
  geom_point() +
  geom_vline(xintercept=0.5, linetype="dashed",
             colour="grey50", linewidth=0.6) +
  scale_colour_manual(values=c(
    "Respiratory"="#d7191c",
    "Cardiovascular"="#2c7bb6",
    "Cardiovascular & Diabetes"="#7b2d8b",
    "Diabetes"="#1a9641"
  )) +
  scale_size_continuous(range=c(3,8), name="PP.H4") +
  scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.2)) +
  labs(title="Open Targets Variant-to-Gene (L2G) Scores",
       subtitle="Point size represents PP.H4 colocalisation probability",
       x="L2G Score", y="Causal Gene (Lead Variant)", colour="Organ") +
  theme_minimal() +
  theme(plot.title=element_text(face="bold", size=11),
        plot.subtitle=element_text(size=9, colour="grey40"),
        legend.position="right",
        axis.text.y=element_text(size=9))

png(file.path(OUT_DIR, "opentargets_L2G_dotplot.png"),
    width=800, height=450, res=120)
print(p4)
dev.off()
cat("Saved: Figure 4 - opentargets_L2G_dotplot.png
")
cat("All figures complete!
")
