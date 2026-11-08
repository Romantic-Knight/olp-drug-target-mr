# Step 1: Extract significant cis-eQTL (P<5e-8) for druggable x eQTLGen genes
# Usage: Rscript extract_eqtl_ivs.R
suppressPackageStartupMessages(library(data.table))

drug_file <- "data/processed/targets/druggable_x_eqtlgen.csv"
eqtlgen_file <- "data/raw/eqtlgen/cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
out_file <- "data/processed/eqtl/eqtlgen_sig_cis_p5e8.csv.gz"

drug <- unique(toupper(trimws(fread(drug_file)$gene)))
cat("druggable x eQTLGen genes:", length(drug), "\n")

cols <- c("SNP", "SNPChr", "SNPPos", "AssessedAllele", "OtherAllele",
          "Zscore", "Pvalue", "GeneSymbol", "GeneChr", "GenePos")
cat("reading eQTLGen (full)...\n")
t0 <- Sys.time()
eq <- fread(eqtlgen_file, select = cols, sep = "\t")
cat("total eQTLGen rows:", nrow(eq), "| read time:", round(difftime(Sys.time(), t0, units="secs"), 1), "s\n")

eq[, GeneSymbol := toupper(trimws(GeneSymbol))]
sig <- eq[GeneSymbol %in% drug & Pvalue < 5e-8]
rm(eq); gc()

cat("significant cis-eQTL rows (P<5e-8):", nrow(sig), "\n")
cat("genes with >=1 significant cis-eQTL:", uniqueN(sig$GeneSymbol), "\n")

setorder(sig, GeneSymbol, Pvalue)
dir.create("data/processed/eqtl", recursive = TRUE, showWarnings = FALSE)
fwrite(sig, out_file, compress = "gzip")
cat("written:", out_file, "\n")

# per-gene IV count distribution
cnt <- sig[, .N, by = GeneSymbol]
cat("\n--- per-gene significant SNP counts ---\n")
cat("quantiles:", paste(round(quantile(cnt$N, c(0, .25, .5, .75, .95, 1))), collapse = " / "), "\n")
cat("genes with 1 SNP:", cnt[N == 1, .N], "| >1:", cnt[N > 1, .N], "\n")
cat("top genes by SNP count:\n")
print(head(cnt[order(-N)], 10))
