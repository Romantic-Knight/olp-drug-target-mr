# SuSiE vs coloc.abf PP.H4 comparison figure.
# figure-style compliant: log-log panels per cohort, verdict palette (no red/green
# binary), direct gene labels, PP.H4 = 0.8 guides, values below 10^-4 plotted at the
# floor as open markers (legend explains), PNG (300 dpi) + PDF (cairo_pdf) both written.
# Usage: Rscript plot_susie_coloc.R
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ggrepel)})

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
save_fig <- function(base, plot, width, height) {
  ggsave(paste0(base, ".png"), plot, width = width, height = height, dpi = 300, bg = "white")
  ggsave(paste0(base, ".pdf"), plot, width = width, height = height, device = cairo_pdf)
  cat("saved", paste0(base, ".png"), "and", paste0(base, ".pdf"), "\n")
}

sus <- fread("results/tables/susie_coloc_results.csv")
abf <- fread("results/tables/coloc_olp_results.csv")[, .(gene, cohort, pp_abf = PP.H4)]
d <- merge(sus, abf, by = c("gene", "cohort"))
d[, verdict := fcase(
  is.na(pp_h4_max), "no_FG_CS",
  pp_abf >= 0.8 & pp_h4_max >= 0.8, "confirmed",
  pp_abf >= 0.8 & pp_h4_max < 0.8,  "reversed",
  pp_abf < 0.8  & pp_h4_max >= 0.8, "rescued",
  default = "neither")]
FLOOR <- 1e-4
d[, x := pmax(pp_abf, FLOOR)]
d[, y := pmax(pp_h4_max, FLOOR, na.rm = TRUE)]
d[, below := pp_abf < FLOOR | pp_h4_max < FLOOR]
d[is.na(pp_h4_max), y := FLOOR]
d[is.na(pp_h4_max), verdict := "no_FG_CS"]
d[verdict == "no_FG_CS", below := FALSE]
d[, cohort := factor(cohort, levels = c("WIDE", "INCLAVO"))]

pal <- c(confirmed = "#0072B2", rescued = "#009E73", reversed = "#E69F00",
         neither = "#9AA5B1", no_FG_CS = "#5A5A5A")

p <- ggplot(d, aes(x = x, y = y)) +
  annotate("rect", xmin = 0.8, xmax = Inf, ymin = 0.8, ymax = Inf, fill = "grey75", alpha = 0.25) +
  geom_vline(xintercept = 0.8, linetype = "dashed", colour = "grey40", linewidth = 0.35) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey40", linewidth = 0.35) +
  annotate("text", x = 0.16, y = 0.93, label = "PP.H4 = 0.8", size = 2.4, colour = "grey25") +
  geom_point(aes(colour = verdict, shape = verdict, alpha = below),
             data = d[verdict != "no_FG_CS"], size = 2.1) +
  geom_point(aes(colour = verdict, shape = verdict), data = d[verdict == "no_FG_CS"], size = 2.6) +
  geom_text_repel(aes(label = gene, colour = verdict),
                  size = 2.2, seed = 42, max.overlaps = 25, min.segment.length = 0.2,
                  segment.size = 0.25, segment.colour = "grey60",
                  box.padding = 0.18, force = 1.2) +
  facet_wrap(~ cohort, ncol = 2) +
  scale_x_log10(limits = c(FLOOR * 0.7, 1.4), breaks = 10^seq(-4, 0), labels = c("1e-4", "1e-3", "0.01", "0.1", "1")) +
  scale_y_log10(limits = c(FLOOR * 0.7, 1.4), breaks = 10^seq(-4, 0), labels = c("1e-4", "1e-3", "0.01", "0.1", "1")) +
  scale_colour_manual(values = pal, name = NULL,
                      labels = c(confirmed = "confirmed (abf & SuSiE \u2265 0.8)",
                                 rescued = "SuSiE-only (\u2265 0.8)",
                                 reversed = "abf-only (overturned)",
                                 neither = "neither",
                                 no_FG_CS = "no FinnGen credible set")) +
  scale_shape_manual(values = c(confirmed = 19, rescued = 19, reversed = 19, neither = 19, no_FG_CS = 0), name = NULL,
                     labels = c(confirmed = "confirmed (abf & SuSiE \u2265 0.8)",
                                rescued = "SuSiE-only (\u2265 0.8)",
                                reversed = "abf-only (overturned)",
                                neither = "neither",
                                no_FG_CS = "no FinnGen credible set")) +
  scale_alpha_manual(values = c(`TRUE` = 0.35, `FALSE` = 1), guide = "none") +
  labs(x = "coloc.abf PP.H4 (single-signal model)",
       y = "coloc.susie PP.H4 (max over credible-set pairs)",
       title = "SuSiE multi-signal coloc overturns CD44 & TDP2 and rescues 6 HLA-region genes (WIDE)",
       subtitle = "Open circles = value below 10\u207b\u2074 (plotted at floor); OLP = oral lichen planus, FinnGen R13") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", legend.text = element_text(size = 7),
        plot.title = element_text(size = 9.5, hjust = 0, face = "bold"),
        plot.subtitle = element_text(size = 7.5, hjust = 0, colour = "grey25"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(size = 8.5))

save_fig("figures/fig_susie_vs_abf_pph4", p, 9, 6.4)
cat("Done.\n")
