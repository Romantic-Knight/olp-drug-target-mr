# Step 2: top cis-eQTL per gene (most significant SNP) as instrument
# Usage: Rscript top_cis_eqtl.R
suppressPackageStartupMessages(library(data.table))

in_file <- "data/processed/eqtl/eqtlgen_sig_cis_p5e8.csv.gz"
out_file <- "data/processed/eqtl/top_cis_eqtl.csv"

sig <- fread(in_file)
cat("significant rows:", nrow(sig), "| genes:", uniqueN(sig$GeneSymbol), "\n")

# top SNP per gene: min Pvalue (tie -> lowest position to be deterministic)
setorder(sig, GeneSymbol, Pvalue, SNPPos)
top <- sig[, .SD[1], by = GeneSymbol]
cat("top cis-eQTL genes:", nrow(top), "\n")

# also count how many genes have a top SNP with rsid (not chr format)
top[, is_rs := grepl("^rs", SNP)]
cat("top SNP rsid format:", sum(top$is_rs), "| chr format:", nrow(top) - sum(top$is_rs), "\n")

# keep essential columns for MR
out <- top[, .(gene = GeneSymbol, snp = SNP, chr = SNPChr, pos = SNPPos,
               effect_allele = AssessedAllele, other_allele = OtherAllele,
               zscore = Zscore, pvalue = Pvalue)]
fwrite(out, out_file)
cat("written:", out_file, "| rows:", nrow(out), "\n")

# quick sanity: distribution of p-values
cat("top SNP p-value range:", format(min(out$pvalue), digits=3), "-", format(max(out$pvalue), digits=3), "\n")
cat("genes with p < 1e-20:", out[pvalue < 1e-20, .N], "\n")
