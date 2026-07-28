#!/usr/bin/env Rscript
# ============================================================================
# Script: 08_coloc_analysis.R
# Author: Teresa Ironman
# Date: 2026
# Description: Runs colocalisation analysis (coloc.abf) for all 10
#              Bonferroni-significant regions between telomere length and
#              fibrosis. Fibrosis treated as case-control (s=0.09, ~1:10
#              case:control ratio per Joof et al. 2026). Telomere treated
#              as quantitative trait. PP.H4 >= 0.8 threshold for
#              colocalisation.
# ============================================================================

library(coloc)
library(writexl)

COLOC_DIR <- "/home/tlti1/DISS/coloc"

regions <- list(
  list(organ="respiratory",     gene="TERT",    N=53167,  s=0.09),
  list(organ="cardiovascular",  gene="ZC3HC1",  N=370381, s=0.09),
  list(organ="cardiovascular",  gene="ATXN2",   N=370381, s=0.09),
  list(organ="diabetes",        gene="PABPC4",  N=356334, s=0.09),
  list(organ="diabetes",        gene="MST1R",   N=356334, s=0.09),
  list(organ="diabetes",        gene="SEC61A2", N=356334, s=0.09),
  list(organ="diabetes",        gene="TRMT1",   N=356334, s=0.09),
  list(organ="diabetes",        gene="ATXN2",   N=356334, s=0.09),
  list(organ="intestinalpanc",  gene="DMC1",    N=419221, s=0.09)
  # Note: HHEX excluded - insufficient shared variants between datasets
)

coloc_results <- data.frame()

for (r in regions) {
  tryCatch({
    fib  <- read.table(file.path(COLOC_DIR, paste0(r$organ, "_", r$gene, ".txt")),
                       header=TRUE, sep="	")
    telo <- read.table(file.path(COLOC_DIR, paste0("telomere_", r$gene, ".txt")),
                       header=TRUE, sep="	")

    common  <- intersect(fib$rsid, telo$rs_id)
    fib_m   <- fib[fib$rsid %in% common, ]
    telo_m  <- telo[telo$rs_id %in% common, ]
    fib_m   <- fib_m[order(fib_m$rsid), ]
    telo_m  <- telo_m[order(telo_m$rs_id), ]

    cat("Region:", r$organ, r$gene, "- Shared variants:", length(common), "
")

    # Fibrosis: case-control trait
    d1 <- list(beta=fib_m$beta, varbeta=fib_m$se^2,
               N=r$N, s=r$s, type="cc", snp=fib_m$rsid)

    # Telomere length: quantitative trait
    d2 <- list(beta=telo_m$beta, varbeta=telo_m$standard_error^2,
               N=462666, type="quant",
               MAF=telo_m$effect_allele_frequency, snp=telo_m$rs_id)

    res <- coloc.abf(d1, d2)

    coloc_results <- rbind(coloc_results, data.frame(
      organ = r$organ, gene = r$gene,
      nsnps = res$summary["nsnps"],
      PP.H0 = res$summary["PP.H0.abf"],
      PP.H1 = res$summary["PP.H1.abf"],
      PP.H2 = res$summary["PP.H2.abf"],
      PP.H3 = res$summary["PP.H3.abf"],
      PP.H4 = res$summary["PP.H4.abf"],
      colocalises = ifelse(res$summary["PP.H4.abf"] >= 0.8, "Yes",
                    ifelse(res$summary["PP.H4.abf"] >= 0.7, "Borderline", "No"))
    ))

  }, error=function(e) {
    cat("Error in", r$organ, r$gene, ":", conditionMessage(e), "
")
  })
}

print(coloc_results)
write_xlsx(coloc_results, file.path(COLOC_DIR, "coloc_results_full.xlsx"))
cat("
Results saved to coloc_results_full.xlsx
")
