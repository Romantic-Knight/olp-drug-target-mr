# SuSiE multi-signal coloc sensitivity analysis (coloc.susie) for OLP drug-target MR
# eQTLGen cis-eQTL (quant, z-based) x FinnGen R13 OLP (cc, beta/varbeta), WIDE + INCLAVO.
# Uses 1000G Phase 3 EUR LD matrices (data/processed/coloc/ld_1000g/<gene>_ld.csv.gz).
#
# Method notes
#  * coloc.susie() normally re-fits via runsusie(), whose unbounded maxit*100 retry loop
#    can hang on HLA-region fits. We instead pre-fit each dataset with susieR::susie_rss
#    (bounded max_iter = 1000) and pass the class-"susie" objects to coloc.susie(), which
#    skips re-fitting (niter/convergence flags are recorded).
#  * Both datasets are aligned to the LD coding (1000G ALT dosage): if the study effect
#    allele equals the LD ALT -> sign +1; equals the LD REF -> sign -1; otherwise dropped.
#  * coloc.susie compares the credible sets (CS) of the two susie fits pair-wise
#    (coloc.bf_bf, priors matched to coloc.abf: p1=1e-4, p2=1e-4, p12=1e-5).
#
# Outputs:
#   results/tables/susie_coloc_results.csv  (one row per gene x cohort)
#   results/tables/susie_coloc_pairs.csv    (per CS-pair PP values)
#   results/tables/susie_coloc_cs.csv       (per-CS detail for both datasets)
# Usage: Rscript run_susie_coloc.R
suppressPackageStartupMessages({library(data.table); library(coloc); library(susieR)})

ROOT <- "D:/AI-Wisp-Science/Drug Target MR"
setwd(ROOT)
ld_dir <- file.path(ROOT, "data/processed/coloc/ld_1000g")
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

eqtl <- fread("data/processed/coloc/coloc_eqtl_windows.csv")
fg_wide <- fread("data/processed/coloc/coloc_fg_wide.csv")
fg_incl <- fread("data/processed/coloc/coloc_fg_inclavo.csv")
genes <- unique(eqtl$gene)

meta <- list(WIDE = list(s = 9526 / (9526 + 490660), N = 9526 + 490660),
             INCLAVO = list(s = 2491 / (2491 + 497695), N = 2491 + 497695))
LMAX <- 15            # max number of susie signals
MAX_ITER <- 1000      # bounded susie_rss iterations (avoids runsusie maxit*100 hang)
PRIORS <- list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)   # match coloc.abf defaults

align_sign <- function(effect_allele, ref, alt) {
  fifelse(effect_allele == alt, 1, fifelse(effect_allele == ref, -1, 0))
}

read_ld <- function(g) {
  snps <- fread(file.path(ld_dir, sprintf("%s_snps.csv", g)))
  ld_dt <- fread(file.path(ld_dir, sprintf("%s_ld.csv.gz", g)))
  rs <- ld_dt[[1]]
  LD <- as.matrix(ld_dt[, -1])
  rownames(LD) <- rs
  colnames(LD) <- rs
  LD <- (LD + t(LD)) / 2
  list(snps = snps, LD = LD)
}

# ---- one pre-fitted susie + CS summary ----
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
  cs_tab <- if (n == 0) NULL else rbindlist(lapply(seq_along(cs), function(i) {
    snps_cs <- cs[[i]]
    pips <- fit$pip[snps_cs]
    top <- names(pips)[which.max(pips)]
    data.table(dataset = label, cs_index = i, n_variants = length(snps_cs),
               top_variant = top, top_pip = round(max(pips), 4),
               pip_min = round(min(pips), 4), pip_max = round(max(pips), 4))
  }))
  list(fit = fit, n_cs = n, niter = fit$niter, converged = fit$converged,
       cs_tab = cs_tab, err = "")
}

run_cohort <- function(g, ld, s1f, fg, meta_one, cohort) {
  ld_snps <- ld$snps; LD <- ld$LD
  eq_g <- eqtl[gene == g][, .(snp, effect_allele, other_allele, zscore, pvalue, NrSamples)]
  fg_g <- fg[, .(snp, ref, alt, beta, sebeta, pval)]
  m <- merge(ld_snps, eq_g, by.x = "rsid", by.y = "snp")
  m <- merge(m, fg_g, by.x = "rsid", by.y = "snp")
  m <- m[!is.na(zscore) & !is.na(beta) & !is.na(sebeta) & pval < 1]
  m[, sign_e := align_sign(effect_allele, ref.x, alt.x)]
  m[, sign_f := align_sign(alt.y, ref.x, alt.x)]
  m <- m[sign_e != 0 & sign_f != 0]
  setorder(m, rsid)

  if (nrow(m) < 10) {
    return(data.table(gene = g, cohort = cohort, n_snps = nrow(m), err = "nsnp<10",
                      n_cs_eqtl_full = s1f$n_cs))
  }
  N_eqtl <- round(median(m$NrSamples))
  d1 <- list(snp = m$rsid, beta = m$sign_e * m$zscore, varbeta = rep(1, nrow(m)),
             MAF = m$maf, N = N_eqtl, type = "quant", LD = LD)
  d2 <- list(snp = m$rsid, beta = m$sign_f * m$beta, varbeta = m$sebeta^2,
             type = "cc", s = meta_one$s, N = meta_one$N, LD = LD)
  s1 <- fit_susie(d1, "eqtl")
  s2 <- fit_susie(d2, cohort)
  err <- paste(s1$err, s2$err, sep = "; ")
  if (s1$err == "" && s2$err == "") {
    res <- tryCatch(do.call(coloc.susie, c(list(s1$fit, s2$fit), PRIORS)),
                    error = function(e) e)
    if (inherits(res, "error")) {
      return(data.table(gene = g, cohort = cohort, n_snps = nrow(m), err = paste0("coloc.susie: ", conditionMessage(res)),
                        n_cs_eqtl = s1$n_cs, n_cs_fg = s2$n_cs, n_cs_eqtl_full = s1f$n_cs,
                        niter_eqtl = s1$niter, conv_eqtl = s1$converged,
                        niter_fg = s2$niter, conv_fg = s2$converged))
    }
    if (is.data.table(res) && all(is.na(res$nsnps))) {
      return(data.table(gene = g, cohort = cohort, n_snps = nrow(m), err = "no CS in one/both datasets",
                        n_cs_eqtl = s1$n_cs, n_cs_fg = s2$n_cs, n_cs_eqtl_full = s1f$n_cs,
                        niter_eqtl = s1$niter, conv_eqtl = s1$converged,
                        niter_fg = s2$niter, conv_fg = s2$converged))
    }
    summ <- res$summary
    best <- summ[which.max(PP.H4.abf)]
    out <- data.table(gene = g, cohort = cohort, n_snps = nrow(m),
                      n_cs_eqtl = s1$n_cs, n_cs_fg = s2$n_cs, n_cs_eqtl_full = s1f$n_cs,
                      n_pairs = nrow(summ),
                      pp_h4_max = round(best$PP.H4.abf, 4),
                      pp_h4_2nd = round(sort(summ$PP.H4.abf, decreasing = TRUE)[2], 4),
                      hit1 = best$hit1, hit2 = best$hit2, idx1 = best$idx1, idx2 = best$idx2,
                      niter_eqtl = s1$niter, conv_eqtl = s1$converged,
                      niter_fg = s2$niter, conv_fg = s2$converged, err = "")
    pairs <- cbind(data.table(gene = g, cohort = cohort), summ)
  } else {
    out <- data.table(gene = g, cohort = cohort, n_snps = nrow(m), err = err,
                      n_cs_eqtl = if (!is.null(s1$n_cs)) s1$n_cs else NA,
                      n_cs_fg = if (!is.null(s2$n_cs)) s2$n_cs else NA,
                      n_cs_eqtl_full = s1f$n_cs)
    pairs <- NULL
  }
  cs <- rbind(s1f$cs_tab, s1$cs_tab, s2$cs_tab)
  if (!is.null(cs)) cs[, `:=`(gene = g, cohort = cohort)]
  list(summary = out, pairs = pairs, cs = cs)
}

out_list <- list(); pair_list <- list(); cs_list <- list()
for (g in genes) {
  cat("== ", g, " ==\n"); flush(stdout())
  ld_file <- file.path(ld_dir, sprintf("%s_ld.csv.gz", g))
  if (!file.exists(ld_file)) { cat("  (no LD)\n"); flush(stdout()); next }
  ld <- read_ld(g)
  eq_g <- eqtl[gene == g][, .(snp, effect_allele, other_allele, zscore, pvalue, NrSamples)]
  mf <- merge(ld$snps, eq_g, by.x = "rsid", by.y = "snp")
  mf <- mf[!is.na(zscore)]
  mf[, sign_e := align_sign(effect_allele, ref, alt)]
  mf <- mf[sign_e != 0]
  d1f <- list(snp = mf$rsid, beta = mf$sign_e * mf$zscore, varbeta = rep(1, nrow(mf)),
              MAF = mf$maf, N = round(median(mf$NrSamples)), type = "quant", LD = ld$LD)
  s1f <- fit_susie(d1f, "eqtl_full")
  cat("   eqtl_full: CS =", s1f$n_cs, " niter =", s1f$niter, " conv =", s1f$converged, "\n"); flush(stdout())
  for (cf in c("WIDE", "INCLAVO")) {
    fg <- if (cf == "WIDE") fg_wide else fg_incl
    r <- tryCatch(run_cohort(g, ld, s1f, fg, meta[[cf]], cf),
                  error = function(e) data.table(gene = g, cohort = cf, err = conditionMessage(e)))
    if (!is.null(r$summary)) out_list[[length(out_list) + 1]] <- r$summary
    else if (is.data.table(r) && "gene" %in% names(r)) out_list[[length(out_list) + 1]] <- r
    if (!is.null(r$pairs)) pair_list[[length(pair_list) + 1]] <- r$pairs
    if (!is.null(r$cs)) cs_list[[length(cs_list) + 1]] <- r$cs
    cat("   ", cf, ": ", if (!is.null(r$summary)) paste("n_snps", r$summary$n_snps,
        "pairs", r$summary$n_pairs, "PP4max", round(r$summary$pp_h4_max, 3)) else paste("err:", r$err), "\n")
    flush(stdout())
  }
}
out <- rbindlist(out_list, fill = TRUE)
fwrite(out, "results/tables/susie_coloc_results.csv")
if (length(pair_list)) fwrite(rbindlist(pair_list, fill = TRUE), "results/tables/susie_coloc_pairs.csv")
if (length(cs_list))   fwrite(rbindlist(cs_list, fill = TRUE), "results/tables/susie_coloc_cs.csv")
cat("\nwritten results/tables/susie_coloc_results.csv (", nrow(out), "rows)\n")
print(out[, .(gene, cohort, n_snps, n_cs_eqtl, n_cs_fg, n_cs_eqtl_full, niter_eqtl, conv_eqtl, pp_h4_max, err)])
