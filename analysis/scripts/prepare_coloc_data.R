# Prepare coloc input: 22 WIDE-Bonferroni genes x cis window (eQTLGen GenePos +/- 1 Mb)
# eQTLGen (GRCh37) cis-window SNPs + FinnGen R13 OLP (WIDE & INCLAVO) matched by rsid.
# Usage: Rscript prepare_coloc_data.R
suppressPackageStartupMessages(library(data.table))

dir.create("data/processed/coloc", recursive = TRUE, showWarnings = FALSE)

# ---- 22 WIDE Bonferroni genes ----
mr <- fread("results/tables/mr_olp_wide_results.csv")
genes <- unique(mr[significant_bonf == TRUE]$gene)
cat("coloc genes (n =", length(genes), "):", paste(genes, collapse = ", "), "\n")

# ---- eQTLGen full (one read), cis window +/- 1 Mb around each gene ----
cat("reading raw eQTLGen...\n")
t0 <- Sys.time()
eqtl <- fread("data/raw/eqtlgen/cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz",
              sep = "\t", header = TRUE,
              select = c("SNP", "SNPChr", "SNPPos", "AssessedAllele", "OtherAllele",
                         "Zscore", "Pvalue", "GeneSymbol", "GeneChr", "GenePos", "NrSamples"),
              colClasses = list(character = c("SNP", "AssessedAllele", "OtherAllele", "GeneSymbol")))
cat("  read", nrow(eqtl), "rows in", round(difftime(Sys.time(), t0, units = "secs"), 1), "s\n")

win <- 1e6
l <- lapply(genes, function(g) {
  sub <- eqtl[GeneSymbol == g & abs(SNPPos - GenePos) <= win]
  sub[, gene := g]
  sub[, .(gene, snp = SNP, chr = SNPChr, pos = SNPPos, effect_allele = AssessedAllele,
          other_allele = OtherAllele, zscore = Zscore, pvalue = Pvalue, NrSamples)]
})
eqwin <- rbindlist(l)
cat("cis-window eQTL SNPs:", nrow(eqwin), "across", length(genes), "genes\n")
cat("  per-gene SNP range:", min(eqwin[, .N, by = gene]$N), "-", max(eqwin[, .N, by = gene]$N), "\n")
fwrite(eqwin, "data/processed/coloc/coloc_eqtl_windows.csv")
cat("written data/processed/coloc/coloc_eqtl_windows.csv\n")

# ---- FinnGen OLP: read full, subset to window rsids ----
all_rs <- unique(eqwin$snp)
cat("unique window rsids:", length(all_rs), "\n")

load_fg <- function(path, tag) {
  cat("reading", basename(path), "...\n")
  t0 <- Sys.time()
  fg <- fread(path, sep = "\t", header = TRUE,
              select = c("rsids", "ref", "alt", "pval", "beta", "sebeta", "af_alt"),
              colClasses = list(character = c("rsids", "ref", "alt")))
  cat("  read", nrow(fg), "rows in", round(difftime(Sys.time(), t0, units = "secs"), 1), "s\n")
  setnames(fg, "rsids", "snp")
  fg <- fg[!is.na(snp) & snp %in% all_rs]
  fg[, snp := trimws(snp)]
  fg <- fg[!duplicated(snp)]
  fg[, outcome := tag]
  fg
}
fg_wide <- load_fg("data/raw/finngen_r13/finngen_R13_K11_LICHEN_PLANUS_WIDE.gz", "WIDE")
fg_incl <- load_fg("data/raw/finngen_r13/finngen_R13_K11_ORAL_LICHEN_PLANUS_INCLAVO.gz", "INCLAVO")
cat("WIDE window rsids matched:", nrow(fg_wide), "| INCLAVO:", nrow(fg_incl), "\n")
fwrite(fg_wide, "data/processed/coloc/coloc_fg_wide.csv")
fwrite(fg_incl, "data/processed/coloc/coloc_fg_inclavo.csv")
cat("written coloc_fg_wide.csv / coloc_fg_inclavo.csv\n")
