# Step-6 docking results figure: AutoDock Vina best affinity for DSigDB candidates.
# GDF11 receptor = 5E4G (human mature domain, 1.5 A), SUOX = 2A9A chain A (chicken ortholog).
# Faceted bar chart; per-panel sorted factors; -7/-6 reference lines. PNG+PDF.
suppressPackageStartupMessages({library(ggplot2)})

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
d <- read.csv("results/tables/docking_results.csv", stringsAsFactors = FALSE)
d$target <- factor(d$target, levels = c("GDF11", "SUOX"))
d$col <- ifelse(d$best_affinity_kcal_mol < -7, "good", ifelse(d$best_affinity_kcal_mol < -6, "moderate", "weak"))
d <- d[order(d$target, d$best_affinity_kcal_mol), ]
# factor levels: per-panel (concatenated GDF11 then SUOX order); each panel only sees its own rows
d$ligand <- factor(d$ligand, levels = unique(d$ligand))

p <- ggplot(d, aes(x = ligand, y = best_affinity_kcal_mol, fill = col)) +
  geom_hline(yintercept = -7, linetype = "dashed", colour = "grey35", linewidth = 0.3) +
  geom_hline(yintercept = -6, linetype = "dashed", colour = "grey70", linewidth = 0.3) +
  geom_col(width = 0.66) +
  geom_text(aes(label = sprintf("%.2f", best_affinity_kcal_mol)), hjust = 1.15, size = 2.3, colour = "white") +
  facet_wrap(~ target, scales = "free_x", ncol = 2,
             labeller = labeller(target = c(GDF11 = "GDF11 (5E4G, human mature domain)",
                                            SUOX = "SUOX (2A9A chain A, chicken ortholog)"))) +
  coord_flip() +
  scale_fill_manual(values = c(good = "#0072B2", moderate = "#E69F00", weak = "#D55E00"),
                    name = NULL, labels = c(good = "< -7 (strong)", moderate = "-6 to -7", weak = "> -6 (weak)")) +
  labs(x = NULL, y = "best AutoDock Vina affinity (kcal/mol)",
       title = "Docking of DSigDB drug candidates to GDF11 and SUOX",
       subtitle = "Grid 30\u00d730\u00d730 \u00c5, exhaustiveness 32, 9 modes; SUOX active site at Mo/SO4 midpoint; GDF11 protein centroid") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", legend.text = element_text(size = 7),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 7.5, colour = "grey25"),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA),
        strip.text = element_text(size = 8.5))

ggsave("figures/fig_docking_results.png", p, width = 9.2, height = 6.4, dpi = 300, bg = "white")
ggsave("figures/fig_docking_results.pdf", p, width = 9.2, height = 6.4, device = cairo_pdf)
cat("saved fig_docking_results.png/pdf\n")
