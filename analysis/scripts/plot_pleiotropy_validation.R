# Pleiotropy validation figures:
#  1) fig_pleiotropy_coloc_heatmap : PP.H4 matrix (5 candidate genes x 9 phenos)
#  2) fig_steiger_direction        : r_exp vs r_out by gene (Steiger directionality)
# figure-style: sentence titles, colour-blind-safe palette, PNG(300dpi)+PDF(cairo_pdf).
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(ggrepel)})

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
save_fig <- function(base, plot, width, height) {
  ggsave(paste0(base, ".png"), plot, width = width, height = height, dpi = 300, bg = "white")
  ggsave(paste0(base, ".pdf"), plot, width = width, height = height, device = cairo_pdf)
  cat("saved", paste0(base, ".png/.pdf"), "\n")
}

PHENO_LAB <- c(
  K11_LICHEN_PLANUS_WIDE = "OLP (WIDE)", K11_ORAL_LICHEN_PLANUS_INCLAVO = "OLP (INCLAVO)",
  K11_COELIAC = "Coeliac disease", T1D_WIDE = "Type 1 diabetes",
  E4_HYTHY_AI_STRICT = "Hypothyroidism (autoimmune)", J10_NASALPOLYP = "Nasal polyps",
  L12_LICHENSCLERATROPH = "Lichen sclerosus", E4_GRAVES_STRICT = "Graves disease",
  RHEUMA_SEROPOS_STRICT = "RA (seropositive)", AUTOIMMUNE = "Any autoimmune disease")

# ---------- 1) coloc heatmap ----------
d <- fread("results/tables/pleiotropy_coloc_results.csv")
cand <- c("GDF11", "SUOX", "EGFL8", "TAPBP", "HLA-DPA1")
d <- d[gene %in% cand]
d[, pheno_lab := factor(PHENO_LAB[pheno], levels = unname(PHENO_LAB))]
d[, gene := factor(gene, levels = rev(cand))]
d[, pph4 := pmin(PP.H4, 1)]
d[, strong := PP.H4 >= 0.8]

p1 <- ggplot(d, aes(x = pheno_lab, y = gene, fill = pph4)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", PP.H4),
                colour = ifelse(pph4 > 0.6, "white", "black")), size = 2.6) +
  scale_fill_gradientn(colours = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
                       values = c(0, 0.25, 0.5, 0.8, 1), limits = c(0, 1), name = "PP.H4") +
  scale_colour_identity() +
  labs(x = NULL, y = NULL,
       title = "Chr12 GDF11/SUOX pleiotropy is a shared causal variant; MHC-locus pleiotropy (EGFL8/TAPBP/HLA-DPA1) is not",
       subtitle = "coloc.abf PP.H4 of eQTLGen cis-eQTL vs FinnGen R13 phenotype (OLP = internal control, reproduces Step-3); threshold PP.H4 = 0.8") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 38, hjust = 1, size = 7.5),
        axis.text.y = element_text(size = 9, face = "italic"),
        panel.grid = element_blank(), legend.position = "bottom",
        plot.title = element_text(size = 9.5, face = "bold"),
        plot.subtitle = element_text(size = 7.5, colour = "grey25"))
save_fig("figures/fig_pleiotropy_coloc_heatmap", p1, 8.2, 3.6)

# ---------- 2) Steiger direction ----------
s <- fread("results/tables/steiger_direction_test.csv")
s[, pheno_lab := PHENO_LAB[pheno]]
s[, gene := factor(gene, levels = c("GDF11", "SUOX", "EGFL8", "TAPBP", "HLA-DPA1"))]
s[, dir_ok := direction == "expression->outcome"]
s[, sig := p_steiger < 0.05]
s[, cls := fifelse(dir_ok & sig, "expression -> outcome (p<0.05)",
               fifelse(!dir_ok & sig, "outcome -> expression (p<0.05)", "not significant"))]
s[, cls := factor(cls, levels = c("expression -> outcome (p<0.05)",
                                  "outcome -> expression (p<0.05)", "not significant"))]

p2 <- ggplot(s, aes(x = r_out, y = r_exp)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_point(aes(colour = cls), size = 2.1, alpha = 0.9) +
  geom_text_repel(aes(label = pheno_lab, colour = cls), size = 2.0, seed = 3,
                  max.overlaps = 25, min.segment.length = 0.1, segment.size = 0.2,
                  box.padding = 0.15, force = 1.5, max.time = 6) +
  facet_wrap(~ gene, ncol = 3, scales = "free") +
  scale_colour_manual(values = c("expression -> outcome (p<0.05)" = "#0072B2",
                                 "outcome -> expression (p<0.05)" = "#D55E00",
                                 "not significant" = "#999999"), name = NULL) +
  labs(x = "Correlation of IV with outcome (r_out)",
       y = "Correlation of IV with eQTL expression (r_exp)",
       title = "Steiger test: candidate IVs tag expression before most outcomes, except EGFL8 with coeliac disease (outcome -> expression)",
       subtitle = "Points above the dashed line support expression -> outcome (cis-eQTL MR assumption); Fisher-z p on r_exp vs r_out") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", strip.text = element_text(face = "italic"),
        plot.title = element_text(size = 9.5, face = "bold"),
        plot.subtitle = element_text(size = 7.5, colour = "grey25"))
save_fig("figures/fig_steiger_direction", p2, 8.6, 4.6)
cat("Done.\n")
