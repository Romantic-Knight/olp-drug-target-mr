# Figure: druggable gene screening visualization (Venn + funnel)
# Panel A: 4-source Venn of druggable gene union (n=9,261)
# Panel B: funnel from druggable union to eQTLGen-intersected exposure set
# Usage: Rscript plot_druggable_sources.R
suppressPackageStartupMessages({
  library(data.table)
  library(VennDiagram)
  library(ggplot2)
  library(grid)
})

outdir <- "figures"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---------------- data ----------------
fin <- fread("data/processed/targets/druggable_genes_final.csv")
eq  <- fread("data/processed/targets/eqtlgen_genes.csv")
inter <- fread("data/processed/targets/druggable_x_eqtlgen.csv")

sets <- list(
  Pharos = fin[in_pharos == TRUE, gene],
  DGIdb  = fin[in_dgidb  == TRUE, gene],
  Finan2017 = fin[in_finan2017 == TRUE, gene],
  OpenTargets = fin[in_opentargets == TRUE, gene]
)
n_union   <- nrow(fin)
n_eq      <- nrow(eq)
n_inter   <- nrow(inter)
pct       <- round(100 * n_inter / n_union, 1)

# ---------------- Panel A: Venn (4 sources) ----------------
# PNG (preview) + PDF (vector delivery)
venn_cols <- c("#0072B2", "#E69F00", "#CC79A7", "#009E73")   # Okabe-Ito (CVD-safe)
venn_grob <- venn.diagram(
  x = sets,
  filename = NULL,
  main = "Druggable genes: four-source union",
  main.cex = 2.0, main.fontfamily = "sans",
  fill = venn_cols, alpha = 0.42,
  col = "white", lwd = 2,
  cex = 1.5, fontfamily = "sans",
  cat.cex = 1.5, cat.fontfamily = "sans", cat.dist = 0.04,
  margin = 0.06,
  ext.text = TRUE, ext.line.lwd = 1,
  na = "remove"
)
venn_png <- file.path(outdir, "fig_druggable_venn.png")
png(venn_png, width = 3400, height = 3000, res = 300, bg = "white")
grid.draw(venn_grob)
dev.off()
venn_pdf <- file.path(outdir, "fig_druggable_venn.pdf")
cairo_pdf(venn_pdf, width = 3400 / 300, height = 3000 / 300, bg = "white")
grid.draw(venn_grob)
dev.off()
cat("Panel A written:", venn_png, file.info(venn_png)$size, "bytes |",
    venn_pdf, file.info(venn_pdf)$size, "bytes\n")

# ---------------- Panel B: funnel (druggable -> eQTLGen exposure) ----------------
funnel <- data.table(
  stage = factor(c("eQTLGen cis-eQTL genes (reference)",
                   "Druggable union (4 sources)",
                   "Druggable ∩ eQTLGen (MR exposure)"),
                 levels = c("eQTLGen cis-eQTL genes (reference)",
                            "Druggable union (4 sources)",
                            "Druggable ∩ eQTLGen (MR exposure)")),
  n = c(n_eq, n_union, n_inter),
  lab = c(format(n_eq, big.mark = ","),
          format(n_union, big.mark = ","),
          sprintf("%s (%.1f%%)", format(n_inter, big.mark = ","), pct))
)
funnel[, is_highlight := stage == "Druggable ∩ eQTLGen (MR exposure)"]

p <- ggplot(funnel, aes(x = n, y = stage)) +
  geom_col(fill = ifelse(funnel$is_highlight, "#0072B2", "#BDBDBD"),
           width = 0.62, color = "white", linewidth = 0.4) +
  geom_text(aes(label = lab), hjust = -0.12, size = 5.2,
            color = ifelse(funnel$is_highlight, "#0072B2", "#404040")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.28)),
                     labels = scales::comma) +
  labs(x = "Number of genes", y = NULL,
       title = "Druggable genes with cis-eQTL support in eQTLGen") +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_line(color = "grey30"),
    axis.ticks.length.x = unit(0.15, "cm"),
    axis.line.x = element_line(color = "grey30"),
    axis.title.x = element_text(margin = margin(t = 8)),
    plot.title = element_text(size = 15, face = "bold", hjust = 0),
    plot.margin = margin(10, 14, 10, 10)
  )

funnel_file <- file.path(outdir, "fig_druggable_funnel.png")
ggsave(funnel_file, plot = p, width = 8.4, height = 3.6, dpi = 300, bg = "white")
funnel_pdf <- file.path(outdir, "fig_druggable_funnel.pdf")
ggsave(funnel_pdf, plot = p, width = 8.4, height = 3.6, device = cairo_pdf)
cat("Panel B written:", funnel_file, file.info(funnel_file)$size, "bytes |",
    funnel_pdf, file.info(funnel_pdf)$size, "bytes\n")
cat("n_union:", n_union, "| n_eq:", n_eq, "| n_inter:", n_inter, "| pct:", pct, "\n")
