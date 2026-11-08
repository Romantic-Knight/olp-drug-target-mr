# PheWAS figure for 5 candidate-gene IVs (GDF11, SUOX, EGFL8, TAPBP, HLA-DPA1).
# FinnGen R13 all-phenotype associations fetched via the official variant API
# (2754 phenotypes per IV). Points = -log10(pval) per phenotype, ordered by
# FinnGen category; colour = direction of the alt-allele effect relative to OLP
# (same / opposite / OLP endpoints themselves); GW threshold 5e-8.
# figure-style: sentence title, no red/green, direct labels via ggrepel,
# PNG (300 dpi) + PDF (cairo_pdf) both written.
# Usage: Rscript plot_phewas_olp.R
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ggrepel)})

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
save_fig <- function(base, plot, width, height) {
  ggsave(paste0(base, ".png"), plot, width = width, height = height, dpi = 300, bg = "white")
  ggsave(paste0(base, ".pdf"), plot, width = width, height = height, device = cairo_pdf)
  cat("saved", paste0(base, ".png"), "and", paste0(base, ".pdf"), "\n")
}

P_GW <- 5e-8
OLP_CODES <- c("K11_LICHEN_PLANUS_WIDE", "K11_ORAL_LICHEN_PLANUS_INCLAVO")
# alt-allele OLP beta (FinnGen alt basis, WIDE) for direction reference
olp_beta <- c(GDF11 = -0.0721463, SUOX = -0.0726783, EGFL8 = 0.19092,
              TAPBP = -0.0958103, `HLA-DPA1` = 0.143372)

d <- fread("results/tables/phewas_olp_all.csv")
d <- d[!is.na(pval) & pval < 1]
d[, mlogp := -log10(pval)]
d[, x := seq_len(.N), by = .(gene)]   # order preserved from file (category_index, phenocode)

d[, is_olp := phenocode %in% OLP_CODES]
d[, dir := fifelse(is_olp, "olp",
              fifelse(sign(beta) == sign(olp_beta[gene]), "same", "opp"))]
d[, sig := pval < P_GW]
d[, plot_class := fifelse(is_olp, "OLP",
                     fifelse(sig & dir == "same", "GW same-dir",
                         fifelse(sig & dir == "opp", "GW opp-dir", "not GW")))]
d[, gene := factor(gene, levels = c("GDF11", "SUOX", "EGFL8", "TAPBP", "HLA-DPA1"))]

# direct labels: OLP endpoints + top-5 non-OLP GW hits per gene
lab <- d[plot_class != "not GW"]
top5 <- d[plot_class %in% c("GW same-dir", "GW opp-dir")][order(gene, pval)][, head(.SD, 5), by = gene]
lab <- unique(rbind(lab[is_olp == TRUE], top5), by = c("gene", "phenocode"))

pal <- c("OLP" = "#000000", "GW same-dir" = "#0072B2", "GW opp-dir" = "#E69F00",
         "not GW" = "#C8CDD2")
shp <- c("OLP" = 18, "GW same-dir" = 16, "GW opp-dir" = 16, "not GW" = 16)

p <- ggplot(d, aes(x = x, y = mlogp)) +
  geom_hline(yintercept = -log10(P_GW), linetype = "dashed", colour = "grey45", linewidth = 0.35) +
  geom_point(aes(colour = plot_class, shape = plot_class), size = 0.55, alpha = 0.75) +
  geom_point(data = d[plot_class == "OLP"], aes(colour = plot_class, shape = plot_class),
             size = 2.2) +
  geom_text_repel(data = lab, aes(label = phenostring, colour = plot_class),
                  size = 2.0, seed = 7, max.overlaps = 20, min.segment.length = 0.15,
                  segment.size = 0.22, segment.colour = "grey55",
                  box.padding = 0.15, force = 1.1, max.time = 8) +
  facet_wrap(~ gene, ncol = 5, scales = "fixed") +
  scale_x_continuous(breaks = c(1, 1000, 2000), limits = c(0, 2800), expand = expansion(mult = c(0.005, 0.01))) +
  scale_y_continuous(limits = c(0, 325), breaks = c(0, 50, 100, 150, 200, 250, 300)) +
  scale_colour_manual(values = pal, name = NULL,
                      labels = c(OLP = "OLP endpoints (WIDE / INCLAVO)",
                                 `GW same-dir` = "GW hit, same dir. as OLP",
                                 `GW opp-dir` = "GW hit, opposite dir.",
                                 `not GW` = "not genome-wide significant")) +
  scale_shape_manual(values = shp, name = NULL,
                     labels = c(OLP = "OLP endpoints (WIDE / INCLAVO)",
                                `GW same-dir` = "GW hit, same dir. as OLP",
                                `GW opp-dir` = "GW hit, opposite dir.",
                                `not GW` = "not genome-wide significant")) +
  labs(x = "FinnGen R13 phenotype (n = 2,754, ordered by disease category)",
       y = "-log10(P)",
       title = "All 5 candidate IVs show genome-wide pleiotropy; chr12 GDF11/SUOX hits are all opposite-direction (autoimmune)",
       subtitle = "Dashed line = P = 5e-8; effect direction of the alt allele relative to OLP; OLP = oral lichen planus, FinnGen R13") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", legend.text = element_text(size = 7),
        plot.title = element_text(size = 9.5, hjust = 0, face = "bold"),
        plot.subtitle = element_text(size = 7.5, hjust = 0, colour = "grey25"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(size = 8.5, face = "italic"),
        axis.text.x = element_text(size = 6.5))

save_fig("figures/fig_phewas_olp_ivs", p, 12, 4.6)
cat("Done.\n")
