# Intersect druggable genes (Pharos + DGIdb + Finan2017 union) with eQTLGen cis-eQTL genes
# Usage: Rscript intersect_eqtlgen_druggable.R
library(data.table)

druggable_file <- "data/processed/targets/druggable_genes_final.csv"
eqtlgen_file   <- "data/raw/eqtlgen/cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
out_intersect  <- "data/processed/targets/druggable_x_eqtlgen.csv"
out_summary    <- "results/tables/druggable_eqtlgen_intersection_summary.csv"

# 1) druggable union
dg <- fread(druggable_file)
dg_genes <- unique(toupper(trimws(dg$gene)))
cat("druggable union genes:", length(dg_genes), "\n")

# 2) eQTLGen unique GeneSymbol (read only the GeneSymbol column)
cat("reading eQTLGen GeneSymbol column...\n")
eq <- fread(eqtlgen_file, select = "GeneSymbol", sep = "\t")
eq_genes <- unique(toupper(trimws(eq[!is.na(GeneSymbol) & GeneSymbol != "", GeneSymbol])))
cat("eQTLGen unique genes:", length(eq_genes), "\n")
rm(eq); gc()

# 3) intersection
inter <- sort(intersect(dg_genes, eq_genes))
cat("druggable x eQTLGen:", length(inter), "\n")

# 4) outputs
fwrite(data.table(gene = inter), out_intersect)
summary_dt <- data.table(
  source      = c("druggable_union", "eqtlgen", "intersection"),
  n_genes     = c(length(dg_genes), length(eq_genes), length(inter)),
  pct_of_drug = c(100, NA, round(100 * length(inter) / length(dg_genes), 2))
)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
fwrite(summary_dt, out_summary)
cat("written:", out_intersect, "|", out_summary, "\n")
