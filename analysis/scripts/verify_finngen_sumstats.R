# Verify FinnGen R13 OLP sumstats: gzip integrity + header
library(data.table)

files <- c("data/raw/finngen_r13/finngen_R13_K11_LICHEN_PLANUS_WIDE.gz",
           "data/raw/finngen_r13/finngen_R13_K11_ORAL_LICHEN_PLANUS_INCLAVO.gz")

for (f in files) {
  cat("===", basename(f), "===\n")
  # gzip integrity: full decompression test
  t0 <- Sys.time()
  con <- gzfile(f, "rb")
  n_lines <- 0
  while (TRUE) {
    block <- readLines(con, n = 200000)
    if (length(block) == 0) break
    n_lines <- n_lines + length(block)
  }
  close(con)
  cat("gzip OK | total lines:", n_lines, "| decompress time:", round(difftime(Sys.time(), t0, units="secs"), 1), "s\n")

  # header (first 3 lines)
  h <- readLines(f, n = 3)
  cat("--- header ---\n")
  cat(h[1], "\n")
  cat(h[2], "\n")
}
