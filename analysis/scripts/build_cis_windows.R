# Build cis windows (±1 Mb) for druggable x eQTLGen genes, merge overlapping (GenomicRanges)
# Usage: Rscript build_cis_windows.R
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
})

genes <- fread("data/processed/targets/druggable_x_eqtlgen.csv")$gene
cat("druggable x eQTLGen genes:", length(genes), "\n")

cat("reading eQTLGen coordinates...\n")
eq <- fread("data/raw/eqtlgen/cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz",
            select = c("GeneSymbol", "GeneChr", "GenePos"), sep = "\t")
eq[, GeneSymbol := toupper(trimws(GeneSymbol))]
eq <- eq[GeneSymbol %in% genes & !is.na(GeneChr) & !is.na(GenePos)]
cat("eQTLGen rows for target genes:", nrow(eq), "\n")

# per-gene probe span, then ±1Mb
span <- eq[, .(start = min(GenePos) - 1e6, end = max(GenePos) + 1e6, n_probes = .N), by = .(chr = GeneChr, gene = GeneSymbol)]
span <- span[start < 1, start := 1]
cat("gene windows:", nrow(span), "\n")

# merge with GenomicRanges
gr <- GRanges(seqnames = as.character(span$chr), ranges = IRanges(start = span$start, end = span$end))
gr <- reduce(gr)
merged <- data.table(chr = as.character(seqnames(gr)), start = start(gr), end = end(gr))
cat("merged cis windows:", nrow(merged), "\n")
cat("total spanned bases:", sum(merged$end - merged$start + 1) / 1e6, "Mb\n")

dir.create("data/processed/targets", recursive = TRUE, showWarnings = FALSE)
fwrite(merged, "data/processed/targets/cis_windows_1mb_merged.tsv")
fwrite(merged, "data/processed/targets/cis_windows_1mb.bed", col.names = FALSE, sep = "\t")
cat("written cis_windows_1mb_merged.tsv + cis_windows_1mb.bed\n")
