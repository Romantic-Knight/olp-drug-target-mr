# SuSiE multi-signal sensitivity re-check of the Step-5b pleiotropy colocalisation.
# The chr12 GDF11/SUOX pleiotropy signals (abf PP.H4 > 0.8 for T1D, hypothyroid,
# nasal polyps, lichen sclerosus, autoimmune) are re-tested with coloc.susie to
# confirm they are not artefacts of single-signal coloc.abf in a multi-signal locus.
#
# Reuses the run_susie_coloc.R machinery: susie_rss pre-fitting (bounded max_iter),
# LD from 1000G Phase 3 EUR, alignment to LD coding, coloc.susie with abf priors.
# FinnGen side = the Step-5b regional files (data/processed/phewas/coloc/fg_*_<gene>.csv);
# OLP WIDE is included as an internal control (must reproduce Step-4 WIDE PP.H4).
#
# Outputs:
#   results/tables/susie_pleiotropy_results.csv   (gene x pheno, abf vs susie PP.H4)
#   results/tables/susie_pleiotropy_pairs.csv     (per CS-pair PP values)
# Usage: Rscript run_susie_pleiotropy.R
suppressPackageStartupMessages({library(data.table); library(coloc); library(susieR)})

ROOT <- "D:/AI-Wisp-Science/Drug Target MR"
setwd(ROOT)
ld_dir <- file.path(ROOT, "data/processed/coloc/ld_1000g")
fg_dir <- "data/processed/phewas/coloc"
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

eqtl <- fread("data/processed/coloc/coloc_eqtl_windows.csv")

# cases / controls (FinnGen R13 manifest; OLP WIDE as internal control)
PHENO_META <- list(
  K11_LICHEN_PLANUS_WIDE = c(cases = 9526, ctrls = 490660),
  K11_COELIAC = c(cases = 5308, ctrls = 476554),
  T1D_WIDE = c(cases = 11197, ctrls = 396409),
  E4_HYTHY_AI_STRICT = c(cases = 56098, ctrls = 358456),
  J10_NASALPOLYP = c(cases = 8793, ctrls = 369050),
  L12_LICHENSCLERATROPH = c(cases = 4181, ctrls = 462830),
  E4_GRAVES_STRICT = c(cases = 4247, ctrls = 495939),
  RHEUMA_SEROPOS_STRICT = c(cases = 6358, ctrls = 488366),
  AUTOIMMUNE = c(cases = 131428, ctrls = 368758))
GENES <- c("GDF11", "SUOX")
PHENOS <- names(PHENO_META)

LMAX <- 15
MAX_ITER <- 1000
PRIORS <- list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)

align_sign <- function(effect_allele, ref, alt) {
  fifelse(effect_allele == alt, 1, fifelse(effect_allele == ref, -1, 0))
}

read_ld <- function(g) {
  snps <- fread(file.path(ld_dir, sprintf("%s_snps.csv", g)))
  ld_dt <- fread(file.path(ld_dir, sprintf("%s_ld.csv.gz", g)))
  rs <- ld_dt[[1]]
  LD <- as.matrix(ld_dt[, -1])
  rownames(LD) <- rs; colnames(LD) <- rs
  LD <- (LD + t(LD)) / 2
  list(snps = snps, LD = LD)
}

fit_susie <- function(d, label) {
  z <- d$beta / sqrt(d$varbeta)
  names(z) <- d$snp
  LDsub <- d$LD[d$snp, d$snp]
  fit <- tryCatch(susieR::susie_rss(z = z, R = LDsub, n = d$N, L = LMAX, max_iter = MAX_ITER),
                  error = function(e) e)
  if (inherits(fit, "error")) return(list(err = conditionMessage(fit)))
  fit <- coloc:::annotate_susie(fit, d$snp, LDsub)
  cs <- fit$sets$cs
  n <- if (is.null(cs)) 0 else length(cs)
  list(fit = fit, n_cs = n, niter = fit$niter, converged = fit$converged, err = "")
}

run_one <- function(g, pheno, ld) {
  if (pheno == "K11_LICHEN_PLANUS_WIDE") {
    fg <- fread("data/processed/coloc/coloc_fg_wide.csv")   # Step-3 internal control
  } else {
    fg <- fread(file.path(fg_dir, sprintf("fg_%s_%s.csv", pheno, g)))
  }
  meta_one <- PHENO_META[[pheno]]
  s <- meta_one["cases"] / (meta_one["cases"] + meta_one["ctrls"])
  N <- meta_one["cases"] + meta_one["ctrls"]
  eq_g <- eqtl[gene == g][, .(snp, effect_allele, other_allele, zscore, pvalue, NrSamples)]
  fg_g <- fg[, .(snp, ref, alt, beta, sebeta, pval)]
  m <- merge(ld$snps, eq_g, by.x = "rsid", by.y = "snp")
  m <- merge(m, fg_g, by.x = "rsid", by.y = "snp")
  m <- m[!is.na(zscore) & !is.na(beta) & !is.na(sebeta) & pval < 1]
  m[, sign_e := align_sign(effect_allele, ref.x, alt.x)]
  m[, sign_f := align_sign(alt.y, ref.x, alt.x)]
  m <- m[sign_e != 0 & sign_f != 0]
  setorder(m, rsid)
  if (nrow(m) < 10) return(list(summary = data.table(gene = g, pheno = pheno, n_snps = nrow(m), err = "nsnp<10")))

  N_eqtl <- round(median(m$NrSamples))
  d1 <- list(snp = m$rsid, beta = m$sign_e * m$zscore, varbeta = rep(1, nrow(m)),
             MAF = m$maf, N = N_eqtl, type = "quant", LD = ld$LD)
  d2 <- list(snp = m$rsid, beta = m$sign_f * m$beta, varbeta = m$sebeta^2,
             type = "cc", s = s, N = N, LD = ld$LD)
  s1 <- fit_susie(d1, "eqtl")
  s2 <- fit_susie(d2, pheno)
  err <- paste(s1$err, s2$err, sep = "; ")
  if (s1$err == "" && s2$err == "") {
    res <- tryCatch(do.call(coloc.susie, c(list(s1$fit, s2$fit), PRIORS)), error = function(e) e)
    if (inherits(res, "error"))
      return(list(summary = data.table(gene = g, pheno = pheno, n_snps = nrow(m), err = paste0("coloc.susie: ", conditionMessage(res)),
                        n_cs_eqtl = s1$n_cs, n_cs_fg = s2$n_cs,
                        conv_eqtl = s1$converged, conv_fg = s2$converged)))
    if (is.data.table(res) && all(is.na(res$nsnps)))
      return(list(summary = data.table(gene = g, pheno = pheno, n_snps = nrow(m), err = "no CS in one/both datasets",
                        n_cs_eqtl = s1$n_cs, n_cs_fg = s2$n_cs,
                        conv_eqtl = s1$converged, conv_fg = s2$converged)))
    summ <- res$summary
    best <- summ[which.max(PP.H4.abf)]
    out <- data.table(gene = g, pheno = pheno, n_snps = nrow(m),
                      n_cs_eqtl = s1$n_cs, n_cs_fg = s2$n_cs, n_pairs = nrow(summ),
                      pp_h4_susie = round(best$PP.H4.abf, 4),
                      pp_h4_2nd = round(sort(summ$PP.H4.abf, decreasing = TRUE)[2], 4),
                      hit1 = best$hit1, hit2 = best$hit2,
                      conv_eqtl = s1$converged, conv_fg = s2$converged, err = "")
    pairs <- cbind(data.table(gene = g, pheno = pheno), summ)
    list(summary = out, pairs = pairs)
  } else {
    out <- data.table(gene = g, pheno = pheno, n_snps = nrow(m), err = err,
                      n_cs_eqtl = s1$n_cs, n_cs_fg = s2$n_cs)
    list(summary = out, pairs = NULL)
  }
}

out_list <- list(); pair_list <- list()
args <- commandArgs(trailingOnly = TRUE)      # optional: pheno1 pheno2 ... to run only these
target_phenos <- if (length(args)) intersect(args, PHENOS) else PHENOS
for (g in GENES) {
  ld <- read_ld(g)
  for (ph in target_phenos) {
    if (file.exists("results/tables/susie_pleiotropy_results.csv")) {
      done <- fread("results/tables/susie_pleiotropy_results.csv")
      if (any(done$gene == g & done$pheno == ph)) { cat("skip", g, ph, "(done)\n"); next }
    }
    cat("==", g, ph, "==\n"); flush(stdout())
    r <- tryCatch(run_one(g, ph, ld),
                  error = function(e) list(summary = data.table(gene = g, pheno = ph, err = conditionMessage(e))))
    out_list[[length(out_list) + 1]] <- r$summary
    if (!is.null(r$pairs)) pair_list[[length(pair_list) + 1]] <- r$pairs
    cat("   ", if (!is.null(r$summary$pp_h4_susie)) sprintf("n_snps %d pairs %d PP4_susie %.4f", r$summary$n_snps, r$summary$n_pairs, r$summary$pp_h4_susie) else paste("err:", r$summary$err), "\n")
  }
}
out <- rbindlist(out_list, fill = TRUE)
if (file.exists("results/tables/susie_pleiotropy_results.csv")) {
  old_out <- fread("results/tables/susie_pleiotropy_results.csv")
  old_out[, `:=`(pp_h4_susie = as.numeric(pp_h4_susie), pp_h4_2nd = as.numeric(pp_h4_2nd),
                 n_snps = as.numeric(n_snps), n_cs_eqtl = as.numeric(n_cs_eqtl),
                 n_cs_fg = as.numeric(n_cs_fg), n_pairs = as.numeric(n_pairs),
                 conv_eqtl = as.logical(conv_eqtl), conv_fg = as.logical(conv_fg))]
  out <- rbindlist(list(old_out, out), fill = TRUE)
}
out <- unique(out, by = c("gene", "pheno"))
out <- out[gene %in% GENES]

# merge abf PP.H4 from Step-5b for direct comparison (handles column-name clash with old rows)
abf <- fread("results/tables/pleiotropy_coloc_results.csv")[gene %in% GENES]
abf <- abf[, .(gene, pheno, pp_h4_abf = PP.H4)]
out <- merge(out, abf, by = c("gene", "pheno"), all.x = TRUE, suffixes = c(".old", ".new"))
out[, pp_h4_abf := fifelse(is.na(pp_h4_abf.old), pp_h4_abf.new, pp_h4_abf.old)]
out[, c("pp_h4_abf.old", "pp_h4_abf.new") := NULL]
setcolorder(out, c("gene", "pheno", "pp_h4_abf", "pp_h4_susie", "n_snps", "n_cs_eqtl", "n_cs_fg", "n_pairs", "pp_h4_2nd", "hit1", "hit2", "conv_eqtl", "conv_fg", "err"))
setorder(out, gene, pheno)
fwrite(out, "results/tables/susie_pleiotropy_results.csv")
if (length(pair_list)) {
  new_pairs <- rbindlist(pair_list, fill = TRUE)
  if (file.exists("results/tables/susie_pleiotropy_pairs.csv")) {
    new_pairs <- rbindlist(list(fread("results/tables/susie_pleiotropy_pairs.csv"), new_pairs), fill = TRUE)
  }
  new_pairs <- unique(new_pairs, by = c("gene", "pheno", "hit1", "hit2"))
  fwrite(new_pairs, "results/tables/susie_pleiotropy_pairs.csv")
}

# ---- internal control: OLP WIDE must reproduce Step-4 WIDE coloc.susie PP.H4 ----
step4 <- fread("results/tables/susie_coloc_results.csv")[cohort == "WIDE" & gene %in% GENES]
chk <- merge(out[pheno == "K11_LICHEN_PLANUS_WIDE", .(gene, pp_h4_susie)],
             step4[, .(gene, pp4_step4 = pp_h4_max)], by = "gene")
chk[, diff := abs(pp_h4_susie - pp4_step4)]
cat("\nOLP internal control vs Step-4 WIDE (max |dPP.H4|):", max(chk$diff), "\n")
stopifnot(max(chk$diff) < 1e-6)

cat("\n===== chr12 pleiotropy: abf vs susie PP.H4 =====\n")
print(out[, .(gene, pheno, pp_h4_abf, pp_h4_susie,
              d = pp_h4_susie - pp_h4_abf,
              n_cs_eqtl, n_cs_fg, conv_fg, err)], nrows = 40)
cat("\nsusie PP.H4>0.8:", sum(out$pp_h4_susie > 0.8, na.rm = TRUE), "of", nrow(out), "\n")
