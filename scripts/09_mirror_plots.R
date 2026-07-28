#!/usr/bin/env Rscript
# ============================================================================
# Script: 09_mirror_plots.R
# Author: Teresa Ironman
# Date: 2026
# Description: Produces mirror plots for all 9 testable colocalisation
#              regions. Fibrosis shown upward, telomere downward (scaled to
#              fibrosis range). Points coloured by LD with lead variant
#              (LDLink LDproxy, GBR population) for colocalised regions.
#              MAF >= 1% filter applied. Zero p-values floored to avoid
#              log10(0) = -Inf.
# ============================================================================

library(ggplot2)

COLOC_DIR <- "/home/tlti1/DISS/coloc"
LD_DIR    <- file.path(COLOC_DIR, "LD")
OUT_DIR   <- file.path(COLOC_DIR, "mirror_plots_LD")
dir.create(OUT_DIR, showWarnings=FALSE)

regions <- list(
  list(organ="respiratory",    gene="TERT",    chr=5,
       fib="respiratory_TERT.txt",       telo="telomere_TERT.txt",
       lead="rs7705526",  ld="rs7705526_LD.txt"),
  list(organ="cardiovascular", gene="ZC3HC1",  chr=7,
       fib="cardiovascular_ZC3HC1.txt",  telo="telomere_ZC3HC1.txt",
       lead="rs11556924", ld="rs11556924_LD.txt"),
  list(organ="cardiovascular", gene="ATXN2",   chr=12,
       fib="cardiovascular_ATXN2.txt",   telo="telomere_ATXN2.txt",
       lead="rs4766578",  ld="rs4766578_LD.txt"),
  list(organ="diabetes",       gene="PABPC4",  chr=1,
       fib="diabetes_PABPC4.txt",        telo="telomere_PABPC4.txt",
       lead="rs3768321",  ld="rs3768321_LD.txt"),
  list(organ="diabetes",       gene="MST1R",   chr=3,
       fib="diabetes_MST1R.txt",         telo="telomere_MST1R.txt",
       lead=NA, ld=NA),
  list(organ="diabetes",       gene="TRMT1",   chr=19,
       fib="diabetes_TRMT1.txt",         telo="telomere_TRMT1.txt",
       lead="rs35601737", ld="rs35601737_LD.txt"),
  list(organ="diabetes",       gene="SEC61A2", chr=10,
       fib="diabetes_SEC61A2.txt",       telo="telomere_SEC61A2.txt",
       lead=NA, ld=NA),
  list(organ="diabetes",       gene="ATXN2",   chr=12,
       fib="diabetes_ATXN2.txt",         telo="telomere_ATXN2.txt",
       lead="rs4766578",  ld="rs4766578_LD.txt"),
  list(organ="intestinalpanc", gene="DMC1",    chr=22,
       fib="intestinalpanc_DMC1.txt",    telo="telomere_DMC1.txt",
       lead=NA, ld=NA)
)

# LD colour scale (LocusZoom style)
ld_cols <- c("Lead"="#7b2d8b","0.8-1.0"="#d7191c","0.6-0.8"="#f46d43",
             "0.4-0.6"="#fee090","0.2-0.4"="#74add1","0.0-0.2"="#4575b4")

assign_ld <- function(rsids, ld_data, lead) {
  sapply(rsids, function(rs) {
    if (!is.na(lead) && rs == lead) return("Lead")
    if (is.null(ld_data)) return("0.0-0.2")
    idx <- which(ld_data$RS_Number == rs)
    if (length(idx) == 0) return("0.0-0.2")
    r2 <- ld_data$R2[idx[1]]
    if (is.na(r2))  return("0.0-0.2")
    if (r2 >= 0.8)  return("0.8-1.0")
    if (r2 >= 0.6)  return("0.6-0.8")
    if (r2 >= 0.4)  return("0.4-0.6")
    if (r2 >= 0.2)  return("0.2-0.4")
    return("0.0-0.2")
  })
}

for (r in regions) {
  tryCatch({
    fib  <- read.table(file.path(COLOC_DIR, r$fib),  header=TRUE, sep="	")
    telo <- read.table(file.path(COLOC_DIR, r$telo), header=TRUE, sep="	")

    # MAF >= 1% filter
    fib  <- fib[fib$eaf >= 0.01 & fib$eaf <= 0.99, ]
    telo <- telo[telo$effect_allele_frequency >= 0.01 &
                 telo$effect_allele_frequency <= 0.99, ]

    # Floor zero p-values (floating point underflow)
    min_p_t <- min(telo$p_value[telo$p_value > 0], na.rm=TRUE)
    telo$p_value[telo$p_value == 0] <- min_p_t
    min_p_f <- min(fib$p[fib$p > 0], na.rm=TRUE)
    fib$p[fib$p == 0] <- min_p_f

    ld_data <- if (!is.na(r$ld))
      read.table(file.path(LD_DIR, r$ld), header=TRUE, sep="	", fill=TRUE)
    else NULL

    # Scale telomere y-axis to fibrosis range
    fib_max  <- max(-log10(fib$p), na.rm=TRUE)
    telo_max <- max(-log10(telo$p_value), na.rm=TRUE)

    fib_plot <- data.frame(
      pos     = fib$position / 1e6,
      y       = -log10(fib$p),
      r2_group = assign_ld(fib$rsid, ld_data, r$lead)
    )
    telo_plot <- data.frame(
      pos     = telo$base_pair_location / 1e6,
      y       = -(-log10(telo$p_value) / telo_max * fib_max),
      r2_group = assign_ld(telo$rs_id, ld_data, r$lead)
    )

    combined <- rbind(fib_plot, telo_plot)
    combined$r2_group <- factor(combined$r2_group,
      levels=c("Lead","0.8-1.0","0.6-0.8","0.4-0.6","0.2-0.4","0.0-0.2"))
    combined <- combined[order(combined$r2_group, decreasing=TRUE), ]

    p <- ggplot(combined, aes(x=pos, y=y, colour=r2_group)) +
      geom_point(alpha=0.7, size=1.5) +
      geom_hline(yintercept=0, colour="black", linewidth=0.6) +
      geom_hline(yintercept=-log10(5e-8), linetype="dashed",
                 colour="grey40", linewidth=0.5) +
      geom_hline(yintercept=-(-log10(5e-8)/telo_max*fib_max),
                 linetype="dashed", colour="grey40", linewidth=0.5) +
      scale_colour_manual(values=ld_cols,
                          name="LD with lead
variant (r²)",
                          labels=c("Lead variant","0.8-1.0","0.6-0.8",
                                   "0.4-0.6","0.2-0.4","<0.2"),
                          drop=FALSE) +
      scale_y_continuous(labels=function(x) abs(round(x, 1))) +
      annotate("text", x=min(combined$pos, na.rm=TRUE),
               y=fib_max*0.9, label="Fibrosis ↑",
               colour="#c0392b", size=3.5, hjust=0, fontface="bold") +
      annotate("text", x=min(combined$pos, na.rm=TRUE),
               y=-fib_max*0.9, label="Telomere ↓ (scaled)",
               colour="#1a5276", size=3.5, hjust=0, fontface="bold") +
      labs(
        title    = paste0(r$gene, " region — ",
                          tools::toTitleCase(r$organ),
                          " fibrosis vs Telomere length"),
        subtitle = paste0("Chr", r$chr,
                          " | MAF ≥1% | Telomere y-axis scaled to fibrosis range",
                          if (!is.na(r$lead))
                            paste0(" | LD with ", r$lead, " (GBR)") else ""),
        x = paste0("Position (Mb, chr", r$chr, ")"),
        y = "-log10(p-value)"
      ) +
      theme_minimal() +
      theme(plot.title    = element_text(face="bold", size=11),
            plot.subtitle = element_text(size=9, colour="grey40"),
            legend.position = "right")

    png(file.path(OUT_DIR,
                  paste0(r$organ, "_", r$gene, "_mirror_LD.png")),
        width=1000, height=600, res=120)
    print(p)
    dev.off()
    cat("Done:", r$organ, r$gene, "
")

  }, error=function(e) {
    cat("Error in", r$organ, r$gene, ":", conditionMessage(e), "
")
  })
}
cat("All mirror plots saved to:", OUT_DIR, "
")
