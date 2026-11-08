# Extract NrSamples (per-gene sample size) for the top cis-eQTL instrument of each gene
# from the raw eQTLGen file, then merge into the MR input tables.
# Usage: Rscript extract_top_eqtl_nsamples.R
suppressPackageStartupMessages(library(data.table))

top <- fread("data/processed/eqtl/top_cis_eqtl.csv")
cat("instruments:", nrow(top), "\n")

raw_path <- "data/raw/eqtlgen/cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
cat("reading raw eQTLGen (selected cols)...\n")
t0 <- Sys.time()
raw <- fread(raw_path, sep = "\t", header = TRUE,
             select = c("SNP", "GeneSymbol", "NrSamples"),
             colClasses = list(character = c("SNP", "GeneSymbol"),
                               integer = "NrSamples"))
cat("  read", nrow(raw), "rows in", round(difftime(Sys.time(), t0, units = "secs"), 1), "s\n")
raw <- raw[!duplicated(paste(GeneSymbol, SNP, sep = "|"))]

m <- merge(top, raw, by.x = c("gene", "snp"), by.y = c("GeneSymbol", "SNP"), all.x = TRUE)
cat("NrSamples matched:", sum(!is.na(m$NrSamples)), "/", nrow(m), "\n")
cat("NrSamples summary:\n")
print(summary(m$NrSamples))

fwrite(m[, .(gene, snp, NrSamples)], "data/processed/eqtl/top_cis_eqtl_nsamples.csv")
cat("written data/processed/eqtl/top_cis_eqtl_nsamples.csv\n")
