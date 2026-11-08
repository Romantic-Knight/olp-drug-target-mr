# Coloc result figure: PP.H4 per gene for WIDE and INCLAVO cohorts
# Strong co-localisation threshold PP.H4 > 0.8 marked.
# figure-style compliant; PNG (preview) + PDF (vector) both written.
# Usage: Rscript plot_coloc_results.R
suppressPackageStartupMessages({library(data.table); library(ggplot2)})

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
save_fig <- function(base, plot, width, height) {
  ggsave(paste0(base, ".png"), plot, width = width, height = height, dpi = 300, bg = "white")
  ggsave(paste0(base, ".pdf"), plot, width = width, height = height, device = cairo_pdf)
  cat("saved", paste0(base, ".png"), "and", paste0(base, ".pdf"), "\n")
}

c <- fread("results/tables/coloc_olp_results.csv")
# merge MR significance for context
mr <- fread("results/tables/mr_olp_wide_results.csv")[significant_bonf == TRUE, .(gene, OR_wide = OR, p_wide = p_MR)]
c <- merge(c, mr, by = "gene", all.x = TRUE)

# order genes by WIDE PP.H4
wide_ord <- c[cohort == "WIDE"][order(PP.H4)]$gene
c[, gene := factor(gene, levels = wide_ord)]
c[, cohort := factor(cohort, levels = c("WIDE", "INCLAVO"))]
c[, strong := PP.H4 > 0.8]

p <- ggplot(c, aes(x = gene, y = PP.H4, fill = strong)) +
  geom_col(width = 0.72, alpha = 0.9) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  annotate("text", x = 21.5, y = 0.83, label = "PP.H4 = 0.8 (strong)", size = 2.8, colour = "grey30", hjust = 1) +
  facet_grid(. ~ cohort) +
  scale_fill_manual(values = c(`FALSE` = "#9AA5B1", `TRUE` = "#D1495B"), name = "PP.H4 > 0.8", labels = c("no", "yes")) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25)) +
  coord_flip() +
  labs(x = NULL, y = "Posterior probability of shared causal variant (PP.H4)",
       title = "Co-localisation of druggable-gene cis-eQTL with OLP (FinnGen R13)") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", plot.title = element_text(size = 9, hjust = 0),
        panel.grid.minor = element_blank())

save_fig("figures/fig_coloc_pph4", p, 9, 6.4)
cat("Done.\n")
