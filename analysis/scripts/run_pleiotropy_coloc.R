# Pleiotropy colocalisation: eQTLGen cis-eQTL (per candidate gene) x FinnGen R13
# top pleiotropic phenotypes (from Step-5 PheWAS), using coloc.abf with the same
# priors and formats as run_coloc.R. OLP WIDE is included as an internal control
# (its PP.H4 must reproduce the Step-3 result).
# Output: results/tables/pleiotropy_coloc_results.csv
suppressPackageStartupMessages({library(data.table); library(coloc)})

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

eqwin <- fread("data/processed/coloc/coloc_eqtl_windows.csv")
olp <- fread("data/processed/coloc/coloc_fg_wide.csv")          # internal control
olp[, pheno := "K11_LICHEN_PLANUS_WIDE"]
olp_meta <- list(K11_LICHEN_PLANUS_WIDE = c(cases = 9526, ctrls = 490660))

# (phenotype) -> cases/controls from finngen_R13_manifest.tsv
pheno_meta <- list(
  K11_COELIAC = c(cases = 5308,  ctrls = 476554),
  T1D_WIDE = c(cases = 11197, ctrls = 396409),
  E4_HYTHY_AI_STRICT = c(cases = 56098, ctrls = 358456),
  J10_NASALPOLYP = c(cases = 8793, ctrls = 369050),
  L12_LICHENSCLERATROPH = c(cases = 4181, ctrls = 462830),
  E4_GRAVES_STRICT = c(cases = 4247, ctrls = 495939),
  RHEUMA_SEROPOS_STRICT = c(cases = 6358, ctrls = 488366),
  AUTOIMMUNE = c(cases = 131428, ctrls = 368758))

genes <- unique(eqwin$gene)
phenos <- names(pheno_meta)

# OLP control: reuse Step-3 file (already merged to eQTL windows)
load_olp <- function(g) olp[olp$pheno == "K11_LICHEN_PLANUS_WIDE"]

load_fg <- function(pheno, g) {
  p <- file.path("data/processed/phewas/coloc", sprintf("fg_%s_%s.csv", pheno, g))
  if (!file.exists(p)) return(NULL)
  fread(p)
}

run_coloc_one <- function(eqtl_sub, fg_sub, meta_one) {
  snps <- intersect(eqtl_sub$snp, fg_sub$snp)
  if (length(snps) < 5) return(data.table(nsnp = length(snps), PP.H0 = NA, PP.H1 = NA, PP.H2 = NA, PP.H3 = NA, PP.H4 = NA, err = "nsnp<5"))
  e <- eqtl_sub[snp %in% snps]
  f <- fg_sub[snp %in% snps]
  setkey(e, snp); setkey(f, snp)
  m <- merge(e, f, by = "snp")
  m[, maf_eqtl := pmin(af_alt, 1 - af_alt)]
  d1 <- list(snp = m$snp, pvalues = pmax(m$pvalue, 1e-300), N = m$NrSamples,
             MAF = m$maf_eqtl, type = "quant")
  d2 <- list(snp = m$snp, beta = m$beta, varbeta = m$sebeta^2,
             type = "cc", s = meta_one["cases"] / (meta_one["cases"] + meta_one["ctrls"]),
             N = meta_one["cases"] + meta_one["ctrls"], MAF = m$af_alt)
  res <- tryCatch(coloc.abf(d1, d2), error = function(e) e)
  if (inherits(res, "error")) return(data.table(nsnp = length(snps), PP.H0 = NA, PP.H1 = NA, PP.H2 = NA, PP.H3 = NA, PP.H4 = NA, err = conditionMessage(res)))
  s <- res$summary
  data.table(nsnp = length(snps),
             PP.H0 = s["PP.H0.abf"], PP.H1 = s["PP.H1.abf"], PP.H2 = s["PP.H2.abf"],
             PP.H3 = s["PP.H3.abf"], PP.H4 = s["PP.H4.abf"], err = "")
}

out <- rbindlist(lapply(genes, function(g) {
  e <- eqwin[gene == g]
  # OLP control (from Step-3 data)
  fc <- load_olp(g)
  rc <- run_coloc_one(e, fc, olp_meta$K11_LICHEN_PLANUS_WIDE)
  rc[, `:=`(gene = g, pheno = "K11_LICHEN_PLANUS_WIDE")]
  rows <- list(rc)
  for (ph in phenos) {
    fg <- load_fg(ph, g)
    if (is.null(fg)) next
    r <- run_coloc_one(e, fg, pheno_meta[[ph]])
    r[, `:=`(gene = g, pheno = ph)]
    rows[[length(rows) + 1]] <- r
  }
  rbindlist(rows)
}))
setcolorder(out, c("gene", "pheno", "nsnp", "PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4", "err"))
out[, PP.H4 := as.numeric(PP.H4)]
fwrite(out, "results/tables/pleiotropy_coloc_results.csv")
cat("written results/tables/pleiotropy_coloc_results.csv (", nrow(out), "rows)\n")

# ---- internal control: OLP WIDE PP.H4 must match Step-3 coloc_olp_results.csv ----
old <- fread("results/tables/coloc_olp_results.csv")[cohort == "WIDE"]
chk <- merge(out[pheno == "K11_LICHEN_PLANUS_WIDE"], old, by = "gene", suffixes = c("_new", "_old"))
chk[, diff := abs(PP.H4_new - PP.H4_old)]
cat("\nOLP internal control (max |dPP.H4|):", max(chk$diff), "\n")
stopifnot(max(chk$diff) < 1e-6)

cat("\n===== pleiotropy coloc (PP.H4>0.8 = shared causal variant) =====\n")
d <- out[order(gene, -PP.H4)]
print(d[, .(gene, pheno, nsnp, PP.H4 = round(PP.H4, 3), PP.H3 = round(PP.H3, 3))], nrows = 60)
cat("\nstrong (PP.H4>0.8):", sum(d$PP.H4 > 0.8, na.rm = TRUE), "of", nrow(d), "\n")
