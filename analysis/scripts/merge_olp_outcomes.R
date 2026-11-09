# Step 3 (v2): merge FinnGen R13 OLP outcomes by rsid (eQTLGen is GRCh37, FinnGen GRCh38 -> use rsid)
# Usage: Rscript merge_olp_outcomes.R
suppressPackageStartupMessages(library(data.table))

top <- fread("data/processed/eqtl/top_cis_eqtl.csv")
cat("instruments:", nrow(top), "\n")
top[, snp := trimws(snp)]

# FinnGen columns: #chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_alt af_alt_cases af_alt_controls
load_olp <- function(path, tag) {
  cat("== ", basename(path), " ==\n", sep = "")
  t0 <- Sys.time()
  olp <- fread(path, sep = "\t", header = TRUE,
               select = c("rsids", "ref", "alt", "pval", "beta", "sebeta", "af_alt"),
               colClasses = list(character = c("rsids", "ref", "alt")))
  cat("  read", nrow(olp), "rows in", round(difftime(Sys.time(), t0, units = "secs"), 1), "s\n")
  setnames(olp, "rsids", "snp")
  olp <- olp[!is.na(snp) & snp != ""]
  olp[, snp := trimws(snp)]
  # keep one row per rsid (first)
  olp <- olp[!duplicated(snp)]
  olp[, outcome := tag]
  olp
}

wide <- load_olp("data/raw/finngen_r13/finngen_R13_K11_LICHEN_PLANUS_WIDE.gz", "WIDE")
incl <- load_olp("data/raw/finngen_r13/finngen_R13_K11_ORAL_LICHEN_PLANUS_INCLAVO.gz", "INCLAVO")

merge_align <- function(olp) {
  m <- merge(top, olp, by = "snp", all.x = TRUE)
  # allele alignment: FinnGen beta is for alt allele; eQTLGen effect allele = effect_allele
  # FinnGen `beta` is the effect of the ALT allele; the exposure effect is per the
  # eQTLGen effect allele (AssessedAllele -> effect_allele). Therefore:
  #   ref == effect_allele ("same") -> effect allele is REF -> beta_out = -beta
  #   alt == effect_allele ("flip") -> effect allele is ALT -> beta_out = +beta
  m[, allele_state := fifelse(ref == effect_allele & alt == other_allele, "same",
                       fifelse(ref == other_allele & alt == effect_allele, "flip", "mismatch"))]
  m[, beta_out := fifelse(allele_state == "same", -beta, beta)]
  m[, se_out := sebeta]
  m[, pval_out := pval]
  m[allele_state == "mismatch", c("beta_out", "se_out", "pval_out", "af_alt") := NA]
  m
}

wide_m <- merge_align(wide)
incl_m <- merge_align(incl)

cat("\n--- WIDE ---\n")
cat("matched (non-NA beta_out):", sum(!is.na(wide_m$beta_out)), "/", nrow(wide_m), "\n")
print(table(wide_m$allele_state))
cat("--- INCLAVO ---\n")
cat("matched (non-NA beta_out):", sum(!is.na(incl_m$beta_out)), "/", nrow(incl_m), "\n")
print(table(incl_m$allele_state))

dir.create("data/processed/eqtl", recursive = TRUE, showWarnings = FALSE)
cols_keep <- c("gene", "snp", "chr", "pos", "effect_allele", "other_allele",
               "zscore", "pvalue", "outcome", "beta_out", "se_out", "pval_out",
               "allele_state", "af_alt", "ref", "alt")
fwrite(wide_m[, ..cols_keep], "data/processed/eqtl/mr_input_wide.csv")
fwrite(incl_m[, ..cols_keep], "data/processed/eqtl/mr_input_inclavo.csv")
cat("written mr_input_wide.csv / mr_input_inclavo.csv\n")

for (tag in c("WIDE", "INCLAVO")) {
  d <- fread(sprintf("data/processed/eqtl/mr_input_%s.csv", tolower(tag)))
  cat(sprintf("%s: analysable genes (non-NA beta_out): %d\n", tag, uniqueN(d[!is.na(beta_out)]$gene)))
}
