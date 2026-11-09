# Step 2: Wald-ratio MR of druggable-gene cis-eQTL (eQTLGen) on OLP (FinnGen R13 WIDE + INCLAVO)
# Single top cis-eQTL IV per gene (strong IV, P<5e-8) -> Wald ratio.
# eQTLGen provides Z only, so beta_exp is recovered as Z * se_exp with
#   se_exp ~= 1/sqrt(2*MAF*(1-MAF)*N)   (MAF from FinnGen af_alt aligned to effect allele; N = NrSamples)
# Test statistic: Z_MR = beta_MR/se_MR = (beta_out/beta_exp)/(se_out/|beta_exp|)
#                       = (beta_out/se_out) * sign(beta_exp) = Z_out * sign(Z_exp)   [exact for strong IV]
# Usage: Rscript mr_wald_ratio.R
suppressPackageStartupMessages(library(data.table))

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# ---- load ----
ns <- fread("data/processed/eqtl/top_cis_eqtl_nsamples.csv")
load_input <- function(tag) {
  d <- fread(sprintf("data/processed/eqtl/mr_input_%s.csv", tolower(tag)))
  d <- merge(d, ns, by = c("gene", "snp"), all.x = TRUE)
  d
}
wide <- load_input("WIDE")
incl <- load_input("INCLAVO")

analyze <- function(d, tag) {
  cat("\n=====", tag, "=====\n")
  d0 <- nrow(d)
  # keep genes with outcome data (mismatch/absent already NA in beta_out)
  d <- d[!is.na(beta_out)]
  cat("outcome-matched:", nrow(d), "/", d0, "\n")
  d <- d[!is.na(NrSamples)]
  cat("with NrSamples:", nrow(d), "\n")

  # effect-allele frequency from FinnGen af_alt (af_alt = ALT allele freq); align to the
  # eQTLGen effect allele: "same" means that allele is REF (p_eff = 1 - af_alt), "flip"
  # means it is ALT (p_eff = af_alt). MAF below is min(p, 1-p), so it is unaffected by
  # this convention; the correction makes the reported effect-allele frequency right.
  d[, p_eff := fifelse(allele_state == "same", 1 - af_alt, af_alt)]
  d[, MAF := pmin(p_eff, 1 - p_eff)]
  # guard against extreme MAF -> unstable beta_exp
  d[MAF < 0.005, MAF := 0.005]
  d[, se_exp := 1 / sqrt(2 * MAF * (1 - MAF) * NrSamples)]
  d[, beta_exp := zscore * se_exp]          # same sign as Z_exp by construction
  stopifnot(all(sign(d$beta_exp) == sign(d$zscore)))

  # Wald ratio
  d[, beta_MR := beta_out / beta_exp]
  d[, se_MR := se_out / abs(beta_exp)]
  d[, Z_MR := (beta_out / se_out) * sign(zscore)]   # == beta_MR/se_MR for strong IV
  d[, p_MR := 2 * pnorm(-abs(Z_MR))]
  d[, OR := exp(beta_MR)]
  d[, OR_lo := exp(beta_MR - 1.96 * se_MR)]
  d[, OR_hi := exp(beta_MR + 1.96 * se_MR)]

  n_tests <- nrow(d)
  d[, p_bonf := pmin(p_MR * n_tests, 1)]
  d[, p_fdr := p.adjust(p_MR, method = "BH")]
  d[, significant_bonf := p_MR < 0.05 / n_tests]
  cat("genes tested:", n_tests, "| Bonferroni threshold:", formatC(0.05 / n_tests, format = "e"), "\n")
  cat("Bonferroni-significant:", sum(d$significant_bonf), "\n")
  cat("FDR<0.05:", sum(d$p_fdr < 0.05), "\n")

  cols <- c("gene", "snp", "chr", "pos", "effect_allele", "other_allele", "allele_state",
            "zscore", "pvalue", "NrSamples", "MAF", "beta_exp", "se_exp",
            "beta_out", "se_out", "pval_out", "beta_MR", "se_MR", "Z_MR", "p_MR",
            "OR", "OR_lo", "OR_hi", "p_bonf", "p_fdr", "significant_bonf")
  setorderv(d, "p_MR")
  fwrite(d[, ..cols], sprintf("results/tables/mr_olp_%s_results.csv", tolower(tag)))
  cat("written results/tables/mr_olp_%s_results.csv\n", tolower(tag))
  d[]
}

r_wide <- analyze(wide, "WIDE")
r_incl <- analyze(incl, "INCLAVO")

# ---- cross-cohort replication: WIDE Bonferroni hits checked in INCLAVO ----
cat("\n===== cross-cohort (WIDE discovery -> INCLAVO replication) =====\n")
hits_wide <- r_wide[significant_bonf == TRUE]
rep <- merge(hits_wide[, .(gene, Z_MR_wide = Z_MR, p_wide = p_MR, beta_MR_wide = beta_MR, OR_wide = OR)],
             r_incl[, .(gene, Z_MR_incl = Z_MR, p_incl = p_MR, beta_MR_incl = beta_MR, OR_incl = OR)],
             by = "gene")
rep[, dir_consistent := sign(Z_MR_wide) == sign(Z_MR_incl)]
rep[, repl_nominal := dir_consistent & p_incl < 0.05]
rep[, repl_bonf := dir_consistent & p_incl < 0.05 / nrow(r_incl)]
setorderv(rep, "p_wide")
cat("WIDE Bonferroni hits:", nrow(rep), "\n")
cat("  direction-consistent in INCLAVO:", sum(rep$dir_consistent), "\n")
cat("  replicated (dir + nominal P<0.05):", sum(rep$repl_nominal), "\n")
cat("  replicated (dir + INCLAVO Bonferroni):", sum(rep$repl_bonf), "\n")
fwrite(rep, "results/tables/mr_olp_hits_crosscohort.csv")

# genes tested in only one cohort (outcome-match differences)
g_wide <- unique(r_wide$gene); g_incl <- unique(r_incl$gene)
cat("genes tested WIDE only:", length(setdiff(g_wide, g_incl)),
    "| INCLAVO only:", length(setdiff(g_incl, g_wide)), "\n")

cat("\nDone.\n")
