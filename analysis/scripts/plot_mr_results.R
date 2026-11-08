# MR results figures: volcano (WIDE) + gene-level Manhattan (WIDE & INCLAVO)
# figure-style compliant: 300 dpi, colour-vision safe, sentence-title panels,
# labelled significant genes via ggrepel, threshold lines, n stated.
# Usage: Rscript plot_mr_results.R
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ggrepel)})

dir.create("figures", recursive = TRUE, showWarnings = FALSE)

# PNG (preview) + PDF (vector delivery) for every figure
save_fig <- function(base, plot, width, height) {
  ggsave(paste0(base, ".png"), plot, width = width, height = height, dpi = 300, bg = "white")
  ggsave(paste0(base, ".pdf"), plot, width = width, height = height, device = cairo_pdf)
  cat("saved", paste0(base, ".png"), "and", paste0(base, ".pdf"), "\n")
}

w <- fread("results/tables/mr_olp_wide_results.csv")
i <- fread("results/tables/mr_olp_inclavo_results.csv")
cc <- fread("results/tables/mr_olp_hits_crosscohort.csv")

nw <- nrow(w); ni <- nrow(i)
bonf_w <- 0.05 / nw; bonf_i <- 0.05 / ni
cat(sprintf("n WIDE=%d (bonf %.2e), n INCLAVO=%d (bonf %.2e)\n", nw, bonf_w, ni, bonf_i))

# chr order for Manhattan (eQTLGen GRCh37 positions; top-SNP position as gene proxy)
chr_order <- c(as.character(1:22), "X", "Y", "MT")
w[, chr := factor(chr, levels = chr_order)]
w <- w[!is.na(chr)]
w[, chr_id := as.integer(chr)]
w[, pos := as.numeric(pos)]
w <- w[order(chr_id, pos)]
# cumulative position per chromosome (numeric to avoid integer overflow)
tmp <- w[, .(maxpos = max(pos)), by = chr_id][order(chr_id)]
tmp[, off := cumsum(shift(as.numeric(maxpos), fill = 0, type = "lag"))]
w <- merge(w, tmp[, .(chr_id, off)], by = "chr_id")
w[, xpos := pos + off]
chr_ticks <- w[, .(x = mean(xpos), chrl = chr[1]), by = chr_id]

# labels: WIDE Bonferroni hits; mark the 4 double-Bonferroni replicated genes
w[, label := ifelse(significant_bonf == TRUE, gene, "")]
w[, repl_class := fifelse(significant_bonf & gene %in% cc[repl_bonf == TRUE]$gene, "repl",
                  fifelse(significant_bonf, "bonf", "ns"))]
cat("WIDE: ns", sum(w$repl_class == "ns"), " bonf", sum(w$repl_class == "bonf"),
    " repl", sum(w$repl_class == "repl"), "\n")

# ---------------- Volcano (WIDE) ----------------
# signed x = MR Z statistic; y = -log10 p; semantic zero at x=0
vol <- ggplot(w, aes(x = Z_MR, y = -log10(p_MR))) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = -log10(bonf_w), linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  annotate("text", x = max(w$Z_MR) * 0.72, y = -log10(bonf_w) + 0.35,
           label = sprintf("Bonferroni α = 0.05/%d", nw), size = 2.8, colour = "grey30") +
  geom_point(aes(colour = repl_class, size = repl_class), alpha = 0.75) +
  geom_text_repel(aes(label = label), size = 2.6, max.overlaps = 40,
                  segment.size = 0.25, segment.colour = "grey55", min.segment.length = 0,
                  box.padding = 0.35, point.padding = 0.25, force = 1.5) +
  scale_colour_manual(values = c(ns = "#B0B0B0", bonf = "#2B6CB0", repl = "#D1495B"),
                      name = "WIDE Bonferroni hit",
                      labels = c(ns = "no", bonf = "yes", repl = "yes + INCLAVO Bonf. replicated"),
                      breaks = c("ns", "bonf", "repl")) +
  scale_size_manual(values = c(ns = 0.8, bonf = 1.6, repl = 2.2), guide = "none") +
  scale_x_continuous(limits = c(min(w$Z_MR) - 0.5, max(w$Z_MR) + 0.5)) +
  labs(x = "MR Wald-ratio Z (per-SD expression effect on OLP)",
       y = expression(-log[10](italic(P))),
       title = "cis-eQTL instruments of 5,040 druggable genes vs oral lichen planus (FinnGen R13 WIDE)") +
  theme_bw(base_size = 9) +
  theme(legend.position.inside = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.key.size = unit(0.35, "cm"), legend.text = element_text(size = 7),
        legend.title = element_text(size = 7.5), plot.title = element_text(size = 8.5, hjust = 0))
save_fig("figures/fig_mr_volcano_wide", vol, 7.2, 5.2)

# ---------------- Manhattan (WIDE) ----------------
manh_w <- ggplot(w, aes(x = xpos, y = -log10(p_MR))) +
  geom_point(aes(colour = repl_class), size = 0.9, alpha = 0.8) +
  geom_hline(yintercept = -log10(bonf_w), linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  annotate("text", x = max(w$xpos) * 0.6, y = -log10(bonf_w) + 0.4,
           label = sprintf("Bonferroni α = 0.05/%d (n = %d genes)", nw, nw), size = 2.8, colour = "grey30") +
  geom_text_repel(aes(label = label), size = 2.4, max.overlaps = 50,
                  segment.size = 0.2, segment.colour = "grey55", min.segment.length = 0,
                  box.padding = 0.3, point.padding = 0.2, force = 1.2) +
  scale_colour_manual(values = c(ns = "#9AA5B1", bonf = "#2B6CB0", repl = "#D1495B"), guide = "none") +
  scale_x_continuous(breaks = chr_ticks$x, labels = chr_ticks$chrl, expand = expansion(mult = c(0.005, 0.01))) +
  labs(x = "Chromosome (top cis-eQTL position, GRCh37)",
       y = expression(-log[10](italic(P))),
       title = "WIDE: MR Wald-ratio P by genomic position") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(size = 8.5, hjust = 0), panel.grid.minor = element_blank())
save_fig("figures/fig_mr_manhattan_wide", manh_w, 9, 4.2)

# ---------------- Manhattan (INCLAVO) ----------------
i2 <- copy(i)
i2[, chr := factor(chr, levels = chr_order)]
i2 <- i2[!is.na(chr)]
i2[, chr_id := as.integer(chr)]
i2[, pos := as.numeric(pos)]
i2 <- i2[order(chr_id, pos)]
tmp2 <- i2[, .(maxpos = max(pos)), by = chr_id][order(chr_id)]
tmp2[, off := cumsum(shift(as.numeric(maxpos), fill = 0, type = "lag"))]
i2 <- merge(i2, tmp2[, .(chr_id, off)], by = "chr_id")
i2[, xpos := pos + off]
chr_ticks2 <- i2[, .(x = mean(xpos), chrl = chr[1]), by = chr_id]
i2[, label := ifelse(significant_bonf == TRUE, gene, "")]
i2[, repl_class := fifelse(significant_bonf, "bonf", "ns")]
manh_i <- ggplot(i2, aes(x = xpos, y = -log10(p_MR))) +
  geom_point(aes(colour = repl_class), size = 0.9, alpha = 0.8) +
  geom_hline(yintercept = -log10(bonf_i), linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  annotate("text", x = max(i2$xpos) * 0.6, y = -log10(bonf_i) + 0.35,
           label = sprintf("Bonferroni α = 0.05/%d (n = %d genes)", ni, ni), size = 2.8, colour = "grey30") +
  geom_text_repel(aes(label = label), size = 2.4, max.overlaps = 20,
                  segment.size = 0.2, segment.colour = "grey55", min.segment.length = 0,
                  box.padding = 0.3, point.padding = 0.2, force = 1.2) +
  scale_colour_manual(values = c(ns = "#9AA5B1", bonf = "#2B6CB0"), guide = "none") +
  scale_x_continuous(breaks = chr_ticks2$x, labels = chr_ticks2$chrl, expand = expansion(mult = c(0.005, 0.01))) +
  labs(x = "Chromosome (top cis-eQTL position, GRCh37)",
       y = expression(-log[10](italic(P))),
       title = "INCLAVO (sensitivity): MR Wald-ratio P by genomic position") +
  theme_bw(base_size = 9) +
  theme(plot.title = element_text(size = 8.5, hjust = 0), panel.grid.minor = element_blank())
save_fig("figures/fig_mr_manhattan_inclavo", manh_i, 9, 4.2)
cat("Done.\n")
