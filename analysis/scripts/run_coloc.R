# Coloc (coloc.abf) of eQTLGen cis-eQTL x FinnGen R13 OLP for 22 WIDE-Bonferroni genes
# WIDE (primary) + INCLAVO (sensitivity). PP.H4 > 0.8 = strong co-localisation.
# Usage: Rscript run_coloc.R
suppressPackageStartupMessages({library(data.table); library(coloc)})

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

eqwin <- fread("data/processed/coloc/coloc_eqtl_windows.csv")
fg_wide <- fread("data/processed/coloc/coloc_fg_wide.csv")
fg_incl <- fread("data/processed/coloc/coloc_fg_inclavo.csv")
genes <- unique(eqwin$gene)
cat("genes:", length(genes), "\n")

# FinnGen sampling fraction (case/control ratio)
meta <- list(WIDE = list(s = 9526 / (9526 + 490660), N = 9526 + 490660),
             INCLAVO = list(s = 2491 / (2491 + 497695), N = 2491 + 497695))

run_coloc_one <- function(eqtl_sub, fg_sub, meta_one) {
  snps <- intersect(eqtl_sub$snp, fg_sub$snp)
  if (length(snps) < 5) return(data.table(nsnp = length(snps), PP.H0 = NA, PP.H1 = NA, PP.H2 = NA, PP.H3 = NA, PP.H4 = NA, err = "nsnp<5"))
  e <- eqtl_sub[snp %in% snps]
  f <- fg_sub[snp %in% snps]
  setkey(e, snp); setkey(f, snp)
  m <- merge(e, f, by = "snp")
  # MAF for the eQTL side is borrowed from FinnGen af_alt (minor-allele frequency is
  # allele-direction independent; same approach as beta_exp recovery in the MR step)
  m[, maf_eqtl := pmin(af_alt, 1 - af_alt)]
  d1 <- list(snp = m$snp, pvalues = pmax(m$pvalue, 1e-300), N = m$NrSamples,
             MAF = m$maf_eqtl, type = "quant")
  d2 <- list(snp = m$snp, beta = m$beta, varbeta = m$sebeta^2,
             type = "cc", s = meta_one$s, N = meta_one$N, MAF = m$af_alt)
  res <- tryCatch(coloc.abf(d1, d2), error = function(e) e)
  if (inherits(res, "error")) return(data.table(nsnp = length(snps), PP.H0 = NA, PP.H1 = NA, PP.H2 = NA, PP.H3 = NA, PP.H4 = NA, err = conditionMessage(res)))
  s <- res$summary
  data.table(nsnp = length(snps),
             PP.H0 = s["PP.H0.abf"], PP.H1 = s["PP.H1.abf"], PP.H2 = s["PP.H2.abf"],
             PP.H3 = s["PP.H3.abf"], PP.H4 = s["PP.H4.abf"], err = "")
}

out <- rbindlist(lapply(genes, function(g) {
  e <- eqwin[gene == g]
  rw <- run_coloc_one(e, fg_wide, meta$WIDE)
  rw[, `:=`(gene = g, cohort = "WIDE")]
  ri <- run_coloc_one(e, fg_incl, meta$INCLAVO)
  ri[, `:=`(gene = g, cohort = "INCLAVO")]
  rbind(rw, ri)
}))
setcolorder(out, c("gene", "cohort", "nsnp", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4", "err"))
out[, PP.H4 := as.numeric(PP.H4)]
fwrite(out, "results/tables/coloc_olp_results.csv")
cat("written results/tables/coloc_olp_results.csv (", nrow(out), "rows)\n")

# ---- summary ----
cat("\n===== coloc summary =====\n")
wide <- out[cohort == "WIDE"][order(-PP.H4)]
cat("WIDE: strong (PP.H4>0.8):", sum(wide$PP.H4 > 0.8, na.rm = TRUE),
    "| moderate (0.5-0.8):", sum(wide$PP.H4 > 0.5 & wide$PP.H4 <= 0.8, na.rm = TRUE), "\n")
print(wide[, .(gene, nsnp, PP.H4 = round(PP.H4, 3))])
incl <- out[cohort == "INCLAVO"]
cat("\nINCLAVO strong:", sum(incl$PP.H4 > 0.8, na.rm = TRUE), "\n")
print(incl[, .(gene, nsnp, PP.H4 = round(PP.H4, 3))])
