# Step-5c review figure: chr12 (GDF11/SUOX) pleiotropy under a multi-signal model.
# x = Step-5b coloc.abf PP.H4 (single-signal), y = Step-5c coloc.susie PP.H4
# (max over credible-set pairs). Same verdict palette / log-log / 0.8 guides /
# 10^-4 floor conventions as fig_susie_vs_abf_pph4; panels = gene.
# PNG (300 dpi, bg = white) + PDF (cairo_pdf). Usage: Rscript plot_susie_pleiotropy.R
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ggrepel)})

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
save_fig <- function(base, plot, width, height) {
  ggsave(paste0(base, ".png"), plot, width = width, height = height, dpi = 300, bg = "white")
  ggsave(paste0(base, ".pdf"), plot, width = width, height = height, device = cairo_pdf)
  cat("saved", paste0(base, ".png"), "and", paste0(base, ".pdf"), "\n")
}

d <- fread("results/tables/susie_pleiotropy_results.csv")

# short, readable phenotype labels (FinnGen codes in parentheses on first use in report)
d[, pheno_short := fcase(
  pheno == "K11_LICHEN_PLANUS_WIDE", "OLP WIDE",
  pheno == "K11_COELIAC",            "Coeliac",
  pheno == "T1D_WIDE",               "T1D",
  pheno == "E4_HYTHY_AI_STRICT",     "Hypothyroid",
  pheno == "J10_NASALPOLYP",         "Nasal polyp",
  pheno == "L12_LICHENSCLERATROPH",  "Lichen sclerosus",
  pheno == "E4_GRAVES_STRICT",       "Graves",
  pheno == "RHEUMA_SEROPOS_STRICT",  "RA (seropos.)",
  pheno == "AUTOIMMUNE",             "Autoimmune",
  default = pheno)]

# verdict identical to Step 4: is.na(susie) -> no_FG_CS (n_cs_fg == 0 in all such rows)
d[, verdict := fcase(
  is.na(pp_h4_susie), "no_FG_CS",
  pp_h4_abf >= 0.8 & pp_h4_susie >= 0.8, "confirmed",
  pp_h4_abf >= 0.8 & pp_h4_susie < 0.8,  "reversed",
  pp_h4_abf < 0.8  & pp_h4_susie >= 0.8, "rescued",
  default = "neither")]
FLOOR <- 1e-4
d[, x := pmax(pp_h4_abf, FLOOR)]
d[, y := pmax(pp_h4_susie, FLOOR, na.rm = TRUE)]
d[, below := pp_h4_abf < FLOOR | pp_h4_susie < FLOOR]
d[is.na(pp_h4_susie), y := FLOOR]
d[verdict == "no_FG_CS", below := FALSE]
d[, gene := factor(gene, levels = c("GDF11", "SUOX"))]

pal <- c(confirmed = "#0072B2", rescued = "#009E73", reversed = "#E69F00",
         neither = "#9AA5B1", no_FG_CS = "#5A5A5A")

p <- ggplot(d, aes(x = x, y = y)) +
  annotate("rect", xmin = 0.8, xmax = Inf, ymin = 0.8, ymax = Inf, fill = "grey75", alpha = 0.25) +
  geom_vline(xintercept = 0.8, linetype = "dashed", colour = "grey40", linewidth = 0.35) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey40", linewidth = 0.35) +
  annotate("text", x = 0.16, y = 0.93, label = "PP.H4 = 0.8", size = 2.4, colour = "grey25") +
  geom_point(aes(colour = verdict, shape = verdict, alpha = below),
             data = d[verdict != "no_FG_CS"], size = 2.3) +
  geom_point(aes(colour = verdict, shape = verdict), data = d[verdict == "no_FG_CS"], size = 2.8) +
  geom_text_repel(aes(label = pheno_short, colour = verdict),
                  size = 2.2, seed = 7, max.overlaps = 20, min.segment.length = 0.15,
                  segment.size = 0.25, segment.colour = "grey60",
                  box.padding = 0.18, force = 1.3, max.time = 8) +
  facet_wrap(~ gene, ncol = 2) +
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
  labs(x = "Step-5b coloc.abf PP.H4 (single-signal model)",
       y = "Step-5c coloc.susie PP.H4 (max over credible-set pairs)",
       title = "SuSiE confirms GDF11's six pleiotropy signals but prunes SUOX hypothyroid & autoimmune; OLP coloc robust in both",
       subtitle = "Open squares = no FinnGen credible set; open circles = value below 10\u207b\u2074 (plotted at floor); OLP = oral lichen planus, FinnGen R13") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", legend.text = element_text(size = 7),
        plot.title = element_text(size = 9.5, hjust = 0, face = "bold"),
        plot.subtitle = element_text(size = 7.5, hjust = 0, colour = "grey25"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(size = 8.5, face = "italic"))

save_fig("figures/fig_susie_pleiotropy_review", p, 8.6, 5.6)
cat("Done.\n")
