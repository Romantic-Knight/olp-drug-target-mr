# Extract unique gene names from DGIdb interactions.tsv (v5.0.11, Dec-2023)
# Usage: Rscript extract_dgidb_genes.R
library(data.table)

infile <- "data/raw/DGIdb/interactions.tsv"
outfile <- "data/processed/targets/dgidb_genes.csv"

# fread skips '#' comment lines automatically; header detected from first non-comment line
d <- fread(infile, sep = "\t", header = TRUE, na.strings = c("", "NA"))
cat("rows:", nrow(d), "cols:", ncol(d), "\n")

# Standard gene symbol column; fall back to gene_claim_name when missing
genes <- unique(toupper(trimws(d[!is.na(gene_name) & gene_name != "", gene_name])))
cat("unique gene_name:", length(genes), "\n")

# gene_name missing rows -> use gene_claim_name
missing <- d[is.na(gene_name) | gene_name == ""]
cat("rows missing gene_name:", nrow(missing), "\n")
if (nrow(missing) > 0) {
  claim <- unique(toupper(trimws(missing[!is.na(gene_claim_name) & gene_claim_name != "", gene_claim_name])))
  genes <- unique(c(genes, claim))
  cat("after adding gene_claim_name fallback:", length(genes), "\n")
}

genes <- sort(genes)
out <- data.table(gene = genes)
fwrite(out, outfile)
cat("written:", outfile, "| total genes:", nrow(out), "\n")
