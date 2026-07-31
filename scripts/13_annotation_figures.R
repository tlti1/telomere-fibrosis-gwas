#!/usr/bin/env Rscript
# ============================================================================
# Script: 13_annotation_figures.R
# Author: Teresa Ironman
# Date: 2026
# Description: Produces annotation figures for Bonferroni-significant variants:
#   1. CADD Phred score dot plot for all 10 lead variants
#   2. g:Profiler pathway enrichment dot plot
#   CADD scores obtained from https://cadd.gs.washington.edu (v1.7, GRCh38)
#   Pathway analysis performed using g:Profiler (GO:MF, FDR < 0.05)
# ============================================================================

library(ggplot2)
OUT_DIR <- "/home/tlti1/DISS/figures"

# ---- CADD scores dot plot --------------------------------------------------
cadd_data <- data.frame(
  gene = c("TERT","ZC3HC1","SH2B3/ATXN2","PABPC4",
           "MST1R","SEC61A2","TRMT1","HHEX","DMC1"),
  organ = c("Respiratory","Cardiovascular","Cardiovascular & Diabetes",
            "Diabetes","Diabetes","Diabetes","Diabetes",
            "Diabetes","Intestinal-Pancreatic"),
  cadd = c(3.307, 5.377, 2.815, 2.945, 2.634, 0.477, 1.959, 0.976, 0.475),
  consequence = c("Intronic","Missense","Intronic","Intronic",
                  "Intronic","Intronic","Intronic","Intronic","Intronic")
)
cadd_data$gene <- factor(cadd_data$gene,
                          levels=cadd_data$gene[order(cadd_data$cadd)])
cadd_data$organ <- factor(cadd_data$organ,
  levels=c("Respiratory","Cardiovascular","Cardiovascular & Diabetes",
           "Diabetes","Intestinal-Pancreatic"))

p_cadd <- ggplot(cadd_data, aes(x=cadd, y=gene,
                                 colour=organ, shape=consequence)) +
  geom_point(size=4) +
  geom_vline(xintercept=20, linetype="dashed",
             colour="red", linewidth=0.6) +
  geom_vline(xintercept=10, linetype="dotted",
             colour="orange", linewidth=0.6) +
  scale_colour_manual(values=c(
    "Respiratory"="#d7191c",
    "Cardiovascular"="#2c7bb6",
    "Cardiovascular & Diabetes"="#7b2d8b",
    "Diabetes"="#1a9641",
    "Intestinal-Pancreatic"="#ff7f00"
  )) +
  scale_shape_manual(values=c("Intronic"=16, "Missense"=17)) +
  scale_x_continuous(limits=c(0, 28), breaks=seq(0, 25, 5)) +
  annotate("text", x=20.3, y=8.5,
           label="CADD ≥ 20
(top 1%)",
           hjust=0, size=3, colour="red") +
  annotate("text", x=10.3, y=8.5,
           label="CADD ≥ 10
(top 10%)",
           hjust=0, size=3, colour="orange") +
  labs(
    title="CADD Phred Scores for Bonferroni-Significant Variants",
    subtitle="Variants scored using CADD v1.7 (GRCh38)",
    x="CADD Phred Score", y="Gene",
    colour="Organ", shape="Consequence"
  ) +
  theme_minimal() +
  theme(
    plot.title=element_text(face="bold", size=11),
    plot.subtitle=element_text(size=9, colour="grey40"),
    legend.position="right",
    axis.text.y=element_text(size=10)
  )

png(file.path(OUT_DIR, "cadd_scores_dotplot.png"),
    width=850, height=500, res=120)
print(p_cadd)
dev.off()
cat("Saved: cadd_scores_dotplot.png
")

# ---- g:Profiler pathway enrichment plot ------------------------------------
# Results from g:Profiler (https://biit.cs.ut.ee/gprofiler/gost)
# Query: 10 Bonferroni-significant genes
# Sources: GO:MF, GO:BP, KEGG, Reactome
# Significance threshold: Benjamini-Hochberg FDR < 0.05
gprofiler_data <- data.frame(
  pathway = c("Template-free RNA
nucleotidyltransferase activity",
              "tRNA (guanine-N2)
dimethyltransferase activity",
              "RNA-directed RNA
polymerase activity",
              "tRNA binding",
              "DNA strand exchange
activity"),
  neg_log10_p = c(1.510, 1.510, 1.510, 1.510, 1.510),
  genes = c("TERT","TRMT1","TERT","TERT, TRMT1","DMC1"),
  n_genes = c(1, 1, 1, 2, 1)
)
gprofiler_data$pathway <- factor(gprofiler_data$pathway,
                                  levels=rev(gprofiler_data$pathway))

p_gp <- ggplot(gprofiler_data,
               aes(x=neg_log10_p, y=pathway, size=n_genes, colour=genes)) +
  geom_point() +
  geom_vline(xintercept=-log10(0.05), linetype="dashed",
             colour="red", linewidth=0.6) +
  scale_size_continuous(range=c(4, 8), name="Genes
in term") +
  scale_colour_manual(values=c(
    "TERT"="#d7191c",
    "TRMT1"="#1a9641",
    "TERT, TRMT1"="#7b2d8b",
    "DMC1"="#ff7f00"
  )) +
  scale_x_continuous(limits=c(0, 2.5)) +
  annotate("text", x=-log10(0.05)+0.05, y=0.6,
           label="FDR = 0.05", hjust=0, size=3, colour="red") +
  labs(
    title="Pathway Enrichment of Bonferroni-Significant Genes",
    subtitle="g:Profiler GO Molecular Function | FDR < 0.05",
    x="-log10(adjusted p-value)",
    y="Pathway", colour="Genes"
  ) +
  theme_minimal() +
  theme(
    plot.title=element_text(face="bold", size=10),
    plot.subtitle=element_text(size=8, colour="grey40"),
    legend.position="right",
    legend.title=element_text(size=9),
    legend.text=element_text(size=8),
    axis.text.y=element_text(size=9),
    plot.margin=margin(10, 10, 10, 10)
  )

png(file.path(OUT_DIR, "gprofiler_pathways.png"),
    width=1000, height=400, res=120)
print(p_gp)
dev.off()
cat("Saved: gprofiler_pathways.png
")
cat("All annotation figures complete!
")
