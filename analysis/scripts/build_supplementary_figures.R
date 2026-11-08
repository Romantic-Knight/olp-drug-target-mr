# build_supplementary_figures.R — S-Fig1..S-Fig4 for the OLP drug-target MR manuscript
# Outputs: figures/supplementary/  (PNG 300 dpi + PDF, 180 mm wide, journal theme)
# 2026-09-13: the supplementary set was reordered to follow the manuscript narrative, so the
# output names no longer follow the code order: the multiplicity panel (built second here)
# is Supplementary Figure 2 and the colocalisation posterior panel (built first) is
# Supplementary Figure 3. Plate content is unchanged.
suppressMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(data.table)
})
source("analysis/scripts/journal_theme.R")

outdir <- "figures/supplementary"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
TBL <- "results/tables"
SUP <- file.path(TBL, "supplementary")
RET <- file.path(SUP, "_retired")   # former tables whose content a figure displays (2026-09-23)

save_supp <- function(base, p, width_mm = 180, height_mm = 120) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  ggplot2::ggsave(file.path(outdir, paste0(base, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(outdir, paste0(base, ".pdf")), p, width = w, height = h, device = grDevices::cairo_pdf)
  cat("saved:", base, "\n")
}

# ---------------------------------------------------------------------------
# S-Fig1 — all 22 WIDE Bonferroni hits: dual-cohort Wald-ratio MR forest
# ---------------------------------------------------------------------------
t1 <- fread(file.path(SUP, "suppl_table_02_wide_hits_mr_replication_coloc.csv"))
# Final curated tier (manuscript): NOT the automated screen tier strings in the table
t1[, class := ifelse(gene == "GDF11", "TIER1",
              ifelse(gene %in% c("SUOX", "EGFL8"), "TIER2",
              ifelse(gene %in% c("TAPBP", "HLA-DPA1"), "TIER3", "other")))]
ord <- t1[order(-log10(MR_P_WIDE)), unique(gene)]
t1[, gene := factor(gene, levels = ord)]

lng1 <- rbind(
  t1[, .(gene, class, OR = MR_OR_WIDE, lo = ORlo_WIDE, hi = ORhi_WIDE, cohort = "WIDE")],
  t1[, .(gene, class, OR = MR_OR_INCLAVO, lo = ORlo_INCLAVO, hi = ORhi_INCLAVO, cohort = "INCLAVO")]
)
p1 <- ggplot(lng1, aes(x = OR, y = gene, colour = class, shape = cohort,
                       xmin = lo, xmax = hi)) +
  geom_vline(xintercept = 1, linetype = "dotted", colour = "grey55", linewidth = 0.35) +
  geom_pointrange(position = position_dodge(width = 0.55), size = 0.28, linewidth = 0.35,
                  fatten = 2.6) +
  scale_x_log10(breaks = c(0.2, 0.5, 1, 2, 5, 10), labels = c("0.2", "0.5", "1", "2", "5", "10")) +
  scale_colour_manual(values = c(TIER1 = "#C8102E", TIER2 = "#005BBB", TIER3 = "#7F7F7F", other = "#B8B8B8"),
                      labels = c(TIER1 = "TIER1", TIER2 = "TIER2", TIER3 = "TIER3", other = "Not tiered"),
                      name = "Final tier") +
  scale_shape_manual(values = c(WIDE = 19, INCLAVO = 21), name = "OLP definition") +
  labs(x = "Odds ratio per 1-SD higher genetically predicted expression (95% CI, log scale)",
       y = NULL,
       title = "All 22 WIDE Bonferroni-significant druggable genes: dual-cohort MR",
       subtitle = "eQTLGen cis-eQTL single-IV Wald ratio -> FinnGen R13 OLP WIDE (9,526) / INCLAVO (2,491); Bonferroni P = 9.9 x 1e-6") +
  theme_journal(8) +
  theme(panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3))
save_supp("supp_fig01_22hits_forest", p1, height_mm = 130)

# ---------------------------------------------------------------------------
# S-Fig3 — coloc.abf full posterior (H0..H4) for 44 gene-definition runs
# ---------------------------------------------------------------------------
t4 <- fread(file.path(RET, "coloc_abf_susie_44runs.csv"))
lng2 <- melt(t4, id.vars = c("gene", "cohort"), measure.vars = c("PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"),
             variable.name = "hypothesis", value.name = "posterior")
lng2[, hypothesis := factor(hypothesis, levels = rev(c("PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4")))]
lng2[, gene := factor(gene, levels = t4[cohort == "WIDE"][order(-PP.H4), gene])]
lng2[, cohort := factor(cohort, levels = c("WIDE", "INCLAVO"))]
hyp_labs <- c("PP.H0" = "H0: no association", "PP.H1" = "H1: eQTL only",
              "PP.H2" = "H2: OLP only", "PP.H3" = "H3: distinct variants",
              "PP.H4" = "H4: shared variant")
p2 <- ggplot(lng2, aes(x = gene, y = posterior, fill = hypothesis)) +
  geom_col(width = 0.82, linewidth = 0.05, colour = "grey25") +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "black", linewidth = 0.35) +
  facet_wrap(~cohort, nrow = 2, scales = "free_y") +
  scale_fill_manual(values = c("PP.H0" = "#D3D3D3", "PP.H1" = "#9FB6D4", "PP.H2" = "#F4C473",
                               "PP.H3" = "#8FAFAF", "PP.H4" = "#C8102E"),
                    labels = hyp_labs, name = "coloc.abf hypothesis") +
  labs(x = NULL, y = "Posterior probability",
       title = "Colocalisation posteriors for all 22 WIDE hits in both OLP definitions",
       subtitle = "coloc.abf full H0-H4 decomposition; dashed line = strong-colocalisation threshold PP.H4 = 0.8; WIDE top, INCLAVO bottom; genes ordered by WIDE PP.H4") +
  theme_journal(8) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 6),
        legend.key.size = unit(0.3, "cm"))
save_supp("supp_fig03_coloc_full_posteriors", p2, height_mm = 130)

# ---------------------------------------------------------------------------
# S-Fig2 — multiplicity sensitivity: Bonferroni vs FDR (5,040 genes)
# ---------------------------------------------------------------------------
t2 <- fread(file.path(SUP, "suppl_table_01_full_mr_screen_5040.csv"))
bonf_thr <- -log10(0.05 / 5040)   # 9.92e-6 -> ~5.003
mk <- function(d, or, p, pfdr, lab) {
  d[, `:=`(logOR = log2(get(or)), mlogp = -log10(get(p)), pf = get(pfdr))]
  d[, cat := ifelse(mlogp > bonf_thr, "Bonferroni", ifelse(pf < 0.05, "FDR only", "not significant"))]
  d[, def := lab]
}
w3 <- mk(copy(t2), "OR_WIDE", "P_WIDE", "PFDR_WIDE", "WIDE (9,526 cases)")
i3 <- mk(copy(t2), "OR_INCLAVO", "P_INCLAVO", "PFDR_INCLAVO", "INCLAVO (2,491 cases)")
lng3 <- rbind(w3, i3, fill = TRUE)
lng3[, def := factor(def, levels = c("WIDE (9,526 cases)", "INCLAVO (2,491 cases)"))]
n_bonf_w <- lng3[def == "WIDE (9,526 cases)"][cat == "Bonferroni", .N]
n_fdr_w  <- lng3[def == "WIDE (9,526 cases)"][cat == "FDR only", .N]
n_bonf_i <- lng3[def == "INCLAVO (2,491 cases)"][cat == "Bonferroni", .N]
n_fdr_i  <- lng3[def == "INCLAVO (2,491 cases)"][cat == "FDR only", .N]
ann <- data.table(def = c("WIDE (9,526 cases)", "INCLAVO (2,491 cases)"),
                  lab = c(sprintf("Bonferroni: %d   FDR-only: %d", n_bonf_w, n_fdr_w),
                          sprintf("Bonferroni: %d   FDR-only: %d", n_bonf_i, n_fdr_i)))
ann[, def := factor(def, levels = levels(lng3$def))]
p3 <- ggplot(lng3, aes(x = logOR, y = mlogp, colour = cat)) +
  geom_point(size = 0.35, alpha = 0.55, shape = 16) +
  geom_hline(yintercept = bonf_thr, linetype = "dashed", colour = "grey25", linewidth = 0.35) +
  geom_text(data = ann, aes(x = -Inf, y = Inf, label = lab), hjust = -0.05, vjust = 1.6,
            size = 2.6, inherit.aes = FALSE, colour = "black", fontface = "bold") +
  facet_wrap(~def, nrow = 1) +
  scale_colour_manual(values = c(Bonferroni = "#C8102E", "FDR only" = "#E69F00",
                                 "not significant" = "#B8B8B8"), name = NULL) +
  labs(x = "log2(OR) per 1-SD higher expression",
       y = expression(-log[10](P[MR])),
       title = "Multiplicity sensitivity across 5,040 druggable genes",
       subtitle = "Red = Bonferroni-significant (P < 9.9 x 1e-6); orange = additionally significant by Benjamini-Hochberg FDR < 0.05; dashed line = Bonferroni threshold") +
  theme_journal(8) +
  theme(legend.position = "bottom")
save_supp("supp_fig02_fdr_sensitivity", p3, height_mm = 85)

# ---------------------------------------------------------------------------
# S-Fig4 — regional plots (FinnGen OLP) around each candidate cis-eQTL IV, LD-coloured
# v2 (QA-fixed): per-gene free x window per row (patchwork 5x2 grid, no facet
# axis-sharing), IV diamond + rsID label clear of the GW line, genome-wide
# 5e-8 dashed line in every panel, single-line r2 legend strip.
# ---------------------------------------------------------------------------
suppressMessages(library(patchwork))
cands <- data.table(gene = c("GDF11", "SUOX", "EGFL8", "TAPBP", "HLA-DPA1"),
                    iv   = c("rs7302200", "rs773107", "rs41268896", "rs3106190", "rs114390018"),
                    chr  = c(12, 12, 6, 6, 6))
eqtl <- fread("data/processed/coloc/coloc_eqtl_windows.csv")
fgw  <- fread("data/processed/coloc/coloc_fg_wide.csv"); fgi <- fread("data/processed/coloc/coloc_fg_inclavo.csv")
fgw[, outcome := "OLP WIDE"]; fgi[, outcome := "OLP INCLAVO"]
fgall <- rbind(fgw[, .(snp, pval, outcome)], fgi[, .(snp, pval, outcome)])

ld_bin <- function(r2) {
  factor(ifelse(is.na(r2), "LD not estimated",
         ifelse(r2 >= 0.8, "r2 >= 0.8", ifelse(r2 >= 0.6, "r2 0.6-0.8",
         ifelse(r2 >= 0.4, "r2 0.4-0.6", "r2 < 0.4")))),
         levels = c("r2 >= 0.8", "r2 0.6-0.8", "r2 0.4-0.6", "r2 < 0.4", "LD not estimated"))
}
lab_pretty <- c("r2 >= 0.8" = "r\u00b2 \u2265 0.8", "r2 0.6-0.8" = "0.6 \u2264 r\u00b2 < 0.8",
                "r2 0.4-0.6" = "0.4 \u2264 r\u00b2 < 0.6", "r2 < 0.4" = "r\u00b2 < 0.4",
                "LD not estimated" = "LD not estimated")
ld_cols <- c("r2 >= 0.8" = "#D7191C", "r2 0.6-0.8" = "#FC8D59", "r2 0.4-0.6" = "#FEE08B",
             "r2 < 0.4" = "#91BFDB", "LD not estimated" = "#B0B0B0")
GW <- -log10(5e-8)

plist <- list()
for (k in seq_len(nrow(cands))) {
  g <- cands[k, gene]; iv <- cands[k, iv]
  e <- eqtl[gene == g]
  ld <- as.matrix(read.csv(gzfile(paste0("data/processed/coloc/ld_1000g/", g, "_ld.csv.gz")),
                          row.names = 1, check.names = FALSE))
  sent <- which(colnames(ld) == iv)
  r2v <- data.table(snp = rownames(ld), r2 = (ld[, sent])^2)
  e <- merge(e, r2v, by = "snp", all.x = TRUE)
  fin <- fgall[snp %in% e$snp]
  fin <- merge(fin, e[, .(snp, pos, r2)], by = "snp", all.x = TRUE)
  fin[, mlogp := -log10(pmax(pval, 1e-320))][, mb := pos / 1e6][, ldcol := ld_bin(r2)]
  fin[, track := factor(outcome, levels = c("OLP WIDE", "OLP INCLAVO"))]
  xr <- range(fin$mb); xpad <- max(diff(xr) * 0.02, 0.008)
  xlim_g <- c(xr[1] - xpad, xr[2] + xpad)
  ylim_g <- c(0, max(max(fin$mlogp) * 1.18, GW * 1.05))
  xbr <- pretty(xlim_g, n = 4); ybr <- pretty(c(0, ylim_g[2]), n = 4)
  for (tt in c("OLP WIDE", "OLP INCLAVO")) {
    d <- fin[track == tt]
    is_left <- tt == "OLP WIDE"
    p <- ggplot(d, aes(x = mb, y = mlogp)) +
      geom_hline(yintercept = GW, linetype = "dashed", colour = "grey45", linewidth = 0.3) +
      geom_point(aes(colour = ldcol), size = 0.55, alpha = 0.85, shape = 16) +
      scale_colour_manual(values = ld_cols, drop = FALSE) +
      scale_x_continuous(limits = xlim_g, breaks = xbr, expand = c(0, 0)) +
      scale_y_continuous(limits = ylim_g, breaks = ybr, expand = c(0, 0)) +
      theme_journal(7) +
      theme(axis.text = element_text(size = 5.4, colour = "black"),
            plot.margin = margin(1, 2, 1, 1, "pt"), legend.position = "none")
    p <- p + labs(y = if (is_left) expression(-log[10](P[OLP])) else NULL, x = NULL)
    ivd <- d[snp == iv]
    if (nrow(ivd) > 0) {
      ly <- min(ivd$mlogp[1] + diff(ylim_g) * 0.09, ylim_g[2] * 0.96)
      p <- p +
        geom_point(data = ivd, shape = 23, size = 2.0, fill = "black", colour = "white", stroke = 0.3) +
        annotate("text", x = ivd$mb[1], y = ly, label = iv, size = 1.9, fontface = "bold")
    }
    plist[[paste(g, tt)]] <- p
  }
}
mk_row2 <- function(k, bottom = FALSE) {
  g <- cands[k, gene]
  pw <- plist[[paste(g, "OLP WIDE")]]; pi2 <- plist[[paste(g, "OLP INCLAVO")]]
  if (bottom) {
    pw <- pw + labs(x = "Position on chromosome (Mb, GRCh37)")
    pi2 <- pi2 + labs(x = "Position on chromosome (Mb, GRCh37)")
  }
  lbl <- ggplot() + xlim(0, 1) + ylim(0, 1) +
    annotate("text", x = 0.5, y = 0.5, label = g, size = 2.7, fontface = "bold") + theme_void()
  lbl + pw + pi2 + plot_layout(widths = c(0.42, 1, 1))
}
hdr_lbl <- ggplot() + xlim(0, 1) + ylim(0, 1) + theme_void()
hdr_w   <- ggplot() + xlim(0, 1) + ylim(0, 1) +
  annotate("text", x = 0.5, y = 0.5, label = "OLP WIDE", size = 2.6, fontface = "bold") + theme_void()
hdr_i   <- ggplot() + xlim(0, 1) + ylim(0, 1) +
  annotate("text", x = 0.5, y = 0.5, label = "OLP INCLAVO", size = 2.6, fontface = "bold") + theme_void()
header  <- hdr_lbl + hdr_w + hdr_i + plot_layout(widths = c(0.42, 1, 1))
lv <- levels(ld_bin(0))
seg <- data.frame(xmin = 1:length(lv) - 0.44, xmax = 1:length(lv) + 0.44, fill = lv)
legend_strip <- ggplot(seg) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = 0.25, ymax = 0.52, fill = fill)) +
  annotate("text", x = 1:length(lv), y = 0.00, label = unname(lab_pretty[lv]), size = 1.9) +
  scale_fill_manual(values = ld_cols, guide = "none") +
  coord_cartesian(xlim = c(0, length(lv) + 0.4), ylim = c(-0.14, 0.6), clip = "off") +
  theme_void()
p4 <- (header + plot_layout(widths = c(0.42, 1, 1))) / mk_row2(1) / mk_row2(2) / mk_row2(3) /
  mk_row2(4) / mk_row2(5, bottom = TRUE) / legend_strip +
  plot_layout(heights = c(0.35, rep(1, 5), 0.45))
save_supp("supp_fig04_regional_plots", p4, height_mm = 165)

cat("All supplementary figures done ->", outdir, "\n")
