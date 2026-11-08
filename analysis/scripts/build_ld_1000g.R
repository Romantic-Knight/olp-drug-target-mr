# Build 1000G Phase 3 (GRCh37, EUR) LD matrices for coloc genes via remote tabix.
# Per gene: cis window -> VCF region -> EUR GT matrix -> LD (Pearson r).
# Output: data/processed/coloc/ld_1000g/<gene>.rds (rsid, pos, ref, alt, LD, maf)
# Usage: Rscript build_ld_1000g.R [gene1 gene2 ...]  (default: all 22)
suppressPackageStartupMessages({library(data.table); library(Rsamtools); library(GenomicRanges)})

dir.create("data/processed/coloc/ld_1000g", recursive = TRUE, showWarnings = FALSE)

# ---- EUR samples ----
panel <- fread("data/raw/1000g/integrated_call_samples_v3.20130502.ALL.panel")
eur <- panel[super_pop == "EUR"]$sample
cat("EUR samples:", length(eur), "\n")

vcf_base <- "https://1000genomes.s3.amazonaws.com/release/20130502/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

eqwin <- fread("data/processed/coloc/coloc_eqtl_windows.csv")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) genes <- args else genes <- unique(eqwin$gene)
cat("genes:", length(genes), "\n")

# VCF sample order per chromosome (extracted locally from the header)
chr_samples <- new.env()
chr_tf <- new.env()

get_chr <- function(chr) {
  if (exists(chr, envir = chr_samples)) return(get(chr, envir = chr_samples))
  s <- fread(sprintf("data/raw/1000g/chr%s_samples.txt", chr), header = FALSE)$V1
  assign(chr, s, envir = chr_samples)
  assign(chr, TabixFile(sprintf(vcf_base, chr)), envir = chr_tf)
  s
}

# parse GT block of one VCF line -> integer vector (0/1/2, NA) for EUR columns
parse_gt <- function(fields, gt_pos, eur_idx) {
  g <- strsplit(fields[10 + eur_idx], ":", fixed = TRUE)
  g <- vapply(g, function(x) x[gt_pos], character(1))
  # g is "0|1" style (3 chars) or "./."
  a <- substr(g, 1, 1); b <- substr(g, 3, 3)
  a[a == "."] <- NA; b[b == "."] <- NA
  suppressWarnings(as.integer(a) + as.integer(b))
}

build_one <- function(gene) {
  sub <- eqwin[gene == gene]
  chr <- as.character(sub$chr[1])
  if (!chr %in% as.character(1:22)) { cat("  skip chr", chr, "\n"); return(NULL) }
  cat("---", gene, "chr", chr, "window snps:", nrow(sub), "---\n")
  s <- get_chr(chr)
  eur_idx <- which(s %in% eur)
  cat("  eur columns:", length(eur_idx), "\n")

  start <- min(sub$pos) - 20000; end <- max(sub$pos) + 20000
  tf <- get(chr, envir = chr_tf)
  x <- scanTabix(tf, param = GRanges(chr, IRanges(start, end)))[[1]]
  cat("  region lines:", length(x), "\n")

  # wanted rsids
  want <- unique(sub$snp)
  out <- vector("list", length(x))
  k <- 0L
  for (line in x) {
    if (substr(line, 1, 1) == "#") next
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    id <- fields[3]
    if (id == "." || !(id %in% want)) next
    alt <- fields[5]
    if (grepl(",", alt, fixed = TRUE)) next            # biallelic only
    fmt <- strsplit(fields[9], ":", fixed = TRUE)[[1]]
    gt_pos <- which(fmt == "GT")[1]
    if (is.na(gt_pos)) next
    gt <- parse_gt(fields, gt_pos, eur_idx)
    k <- k + 1L
    out[[k]] <- list(rsid = id, pos = as.numeric(fields[2]), ref = fields[4], alt = alt, gt = gt)
  }
  out <- out[seq_len(k)]
  cat("  matched biallelic variants:", k, "\n")
  if (k < 10) { cat("  too few; skipped\n"); return(NULL) }

  rsids <- vapply(out, `[[`, character(1), "rsid")
  pos <- vapply(out, `[[`, numeric(1), "pos")
  ref <- vapply(out, `[[`, character(1), "ref")
  alt <- vapply(out, `[[`, character(1), "alt")
  GT <- do.call(rbind, lapply(out, `[[`, "gt"))           # k x n_eur
  maf <- rowMeans(GT, na.rm = TRUE) / 2
  keep <- !is.na(maf) & maf > 0.01 & maf < 0.99
  cat("  after MAF>0.01 filter:", sum(keep), "\n")
  if (sum(keep) < 10) { cat("  too few after filter; skipped\n"); return(NULL) }
  GT <- GT[keep, , drop = FALSE]; rsids <- rsids[keep]; pos <- pos[keep]
  ref <- ref[keep]; alt <- alt[keep]; maf <- maf[keep]

  LD <- cor(t(GT), use = "pairwise.complete.obs")
  rownames(LD) <- colnames(LD) <- rsids
  # sanity: diagonal 1, symmetry
  cat("  LD dim:", paste(dim(LD), collapse = "x"), "| diag ok:", all(abs(diag(LD) - 1) < 1e-8), "\n")

  saveRDS(list(rsid = rsids, pos = pos, ref = ref, alt = alt, maf = maf, LD = LD,
               n_eur = nrow(panel[super_pop == "EUR"])),
          sprintf("data/processed/coloc/ld_1000g/%s.rds", gene))
  cat("  saved ld_1000g/", gene, ".rds (", k, "->", nrow(LD), ")\n", sep = "")
  NULL
}

for (g in genes) {
  tryCatch(build_one(g), error = function(e) cat("  ERROR on", g, ":", conditionMessage(e), "\n"))
}
cat("Done.\n")
