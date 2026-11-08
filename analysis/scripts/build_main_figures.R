# build_main_figures.R —— 按分析逻辑顺序的 7 张组合多面板主图
# 输出 figures/main_figures/fig01..fig07 (png 300dpi + pdf 矢量)
# Fig2-6 复制已定稿 journal 版；Fig1/Fig7 新组合渲染。
suppressMessages({library(ggplot2); library(data.table); library(ggrepel); library(VennDiagram); library(grid); library(patchwork); library(cowplot)})
source("analysis/scripts/journal_theme.R")

MDIR <- "figures/main_figures"; dir.create(MDIR, recursive = TRUE, showWarnings = FALSE)

save_main <- function(base, plot, width_mm, height_mm) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  ggplot2::ggsave(file.path(MDIR, paste0(base, ".png")), plot, width = w, height = h, dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(MDIR, paste0(base, ".pdf")), plot, width = w, height = h, device = grDevices::cairo_pdf)
  cat("saved", base, "\n")
}

# ============================================================
# 子图渲染函数
# ============================================================

# Fig1A: 研究设计漏斗（与 journal_fig1 同几何，去标题改用面板标签）
render_design <- function() {
  xs <- 1:5
  widths  <- c(0.92, 0.80, 0.70, 0.60, 0.52)
  heights <- c(0.95, 0.80, 0.66, 0.54, 0.46)
  yc <- 2.05
  cols <- c("#005BBB", "#1B9E77", "#2B6CB0", "#D1495B", "#C8102E")
  fw <- data.frame(x = xs, xmin = xs - widths/2, xmax = xs + widths/2,
                   ymin = yc - heights/2, ymax = yc + heights/2, fill = cols)
  step_lab <- c("Druggable\ngenome", "eQTLGen\ncis-eQTL\noverlap", "MR screen\n(2 OLP\ndefinitions)",
                "Colocalisation\n(abf + SuSiE)", "Tiered\ncandidates")
  steps <- data.frame(x = xs, step = step_lab,
                      n = c(9261, 5801, 5040, 44, 5),
                      sub = c("9,261 druggable\ngenes (4 sources)", "5,801 genes\nw/ cis-eQTL",
                              "5,040 genes\ntested (Bonf 9.9e-6)", "44 abf + 44\nSuSiE runs",
                              "5 candidates\n(1xT1,2xT2,2xT3)"))
  steps <- cbind(steps, fw[, c("ymin", "ymax")])
  cohorts <- data.frame(x = c(2.5, 3.7), y = 3.98, cohort = c("WIDE 9,526 cases", "INCLAVO 2,491 cases"))
  ann <- data.frame(x = c(1, 2.4, 3.6, 4.75), y = 1.02,
                    txt = c("FinnGen R13\n(GRCh38)", "eQTLGen\n31,684", "1000G LD\n(EUR, n=503)",
                            "pQTL evidence\nINTERVAL/Knight/\ndeCODE 35,652"))
  ggplot() +
    geom_rect(data = fw, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = fw$fill, alpha = 0.30, colour = "grey30", linewidth = 0.35) +
    geom_segment(data = data.frame(x = xs[-5] + widths[-5]/2, xend = xs[-1] - widths[-1]/2, y = yc, yend = yc),
                 aes(x = x, xend = xend, y = y, yend = yend), colour = "grey55", linetype = "dashed", linewidth = 0.28) +
    geom_text(data = steps, aes(x = x, y = ymax + 0.55, label = step), hjust = 0.5, size = 2.6, fontface = "bold", lineheight = 0.85) +
    geom_text(data = steps, aes(x = x, y = yc, label = n), hjust = 0.5, vjust = 0.5, size = 3.0, fontface = "bold", colour = "black") +
    geom_text(data = steps, aes(x = x, y = ymin - 0.05, label = sub), hjust = 0.5, vjust = 1, size = 2.0, colour = "grey25", lineheight = 0.82) +
    geom_segment(data = cohorts, aes(x = x, xend = x, y = 3.98 - 0.13, yend = 2.55), linetype = "dashed", colour = "grey55", linewidth = 0.26) +
    geom_text(data = cohorts, aes(x = x, y = y, label = cohort), hjust = 0.5, size = 2.3, fontface = "bold", colour = "grey15") +
    geom_text(data = ann, aes(x = x, y = y, label = txt), hjust = 0.5, vjust = 1, size = 2.0, colour = "grey30", lineheight = 0.8) +
    coord_cartesian(xlim = c(0.28, 5.72), ylim = c(0.55, 4.32)) +
    labs(title = "Study design: druggable-genome MR of oral lichen planus") +
    theme_journal(8) +
    theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), plot.margin = margin(6, 6, 6, 6))
}

# Fig1B: druggable Venn (cowplot grob)
render_venn <- function() {
  fin <- fread("data/processed/targets/druggable_genes_final.csv")
  sets <- list(Pharos = fin[in_pharos == TRUE, gene], DGIdb = fin[in_dgidb == TRUE, gene],
               Finan2017 = fin[in_finan2017 == TRUE, gene], OpenTargets = fin[in_opentargets == TRUE, gene])
  venn_cols <- c("#0072B2", "#E69F00", "#CC79A7", "#009E73")
  grob <- venn.diagram(x = sets, filename = NULL, main = NULL,
                       fill = venn_cols, alpha = 0.42, col = "white", lwd = 1.1,
                       cex = 0.62, fontfamily = "sans", cat.cex = 0.68, cat.fontfamily = "sans",
                       cat.dist = 0.075, cat.pos = 0, margin = 0.05, ext.text = TRUE, ext.line.lwd = 0.6, na = "remove")
  cowplot::ggdraw() +
    cowplot::draw_grob(grob, x = 0, y = -0.04, width = 1, height = 0.86) +
    cowplot::draw_label("Druggable genes: four-source union", x = 0.5, y = 0.99, vjust = 1,
                        fontface = "bold", size = 8.5, hjust = 0.5) +
    cowplot::draw_label(sprintf("Pharos / DGIdb / Finan 2017 / OpenTargets (union n = %s)", format(nrow(fin), big.mark = ",")),
                        x = 0.5, y = 0.915, vjust = 1, size = 6.8, colour = "grey35", hjust = 0.5)
}

# Fig1C: druggable funnel
render_funnel <- function() {
  fin <- fread("data/processed/targets/druggable_genes_final.csv")
  eq  <- fread("data/processed/targets/eqtlgen_genes.csv")
  inter <- fread("data/processed/targets/druggable_x_eqtlgen.csv")
  n_union <- nrow(fin); n_eq <- nrow(eq); n_inter <- nrow(inter)
  pct <- round(100 * n_inter / n_union, 1)
  funnel <- data.table(
    stage = factor(c("eQTLGen cis-eQTL genes (reference)", "Druggable union (4 sources)",
                     "Druggable \u2229 eQTLGen (MR exposure)"),
                   levels = c("eQTLGen cis-eQTL genes (reference)", "Druggable union (4 sources)",
                              "Druggable \u2229 eQTLGen (MR exposure)")),
    n = c(n_eq, n_union, n_inter),
    lab = c(format(n_eq, big.mark = ","), format(n_union, big.mark = ","),
            sprintf("%s (%.1f%%)", format(n_inter, big.mark = ","), pct)))
  funnel[, is_highlight := stage == "Druggable \u2229 eQTLGen (MR exposure)"]
  ggplot(funnel, aes(x = n, y = stage)) +
    geom_col(fill = ifelse(funnel$is_highlight, "#0072B2", "#BDBDBD"), width = 0.62, color = "white", linewidth = 0.4) +
    geom_text(aes(label = lab), hjust = -0.08, size = 2.8, color = ifelse(funnel$is_highlight, "#0072B2", "#404040")) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.34)), labels = scales::comma) +
    labs(x = "Number of genes", y = NULL,
         title = "Druggable genes with cis-eQTL support in eQTLGen",
         subtitle = sprintf("%s of %s druggable genes (%.1f%%) have an eQTLGen cis-eQTL instrument",
                            format(n_inter, big.mark=","), format(n_union, big.mark=","), pct)) +
    theme_journal(8)
}

# Fig7A: SuSiE pleiotropy (GDF11/SUOX 两 facet 合一面板)
render_susie_pleiotropy <- function() {
  d <- fread("results/tables/susie_pleiotropy_results.csv")
  d[, pheno_short := fcase(
    pheno == "K11_LICHEN_PLANUS_WIDE", "OLP WIDE", pheno == "K11_COELIAC", "Coeliac",
    pheno == "T1D_WIDE", "T1D", pheno == "E4_HYTHY_AI_STRICT", "Hypothyroid",
    pheno == "J10_NASALPOLYP", "Nasal polyp", pheno == "L12_LICHENSCLERATROPH", "Lichen sclerosus",
    pheno == "E4_GRAVES_STRICT", "Graves", pheno == "RHEUMA_SEROPOS_STRICT", "RA (seropos.)",
    pheno == "AUTOIMMUNE", "Autoimmune", default = pheno)]
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
  ggplot(d, aes(x = x, y = y)) +
    annotate("rect", xmin = 0.8, xmax = Inf, ymin = 0.8, ymax = Inf, fill = "grey75", alpha = 0.22) +
    geom_vline(xintercept = 0.8, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    geom_point(aes(colour = verdict, shape = verdict, alpha = below),
               data = d[verdict != "no_FG_CS"], size = 2.0) +
    geom_point(aes(colour = verdict, shape = verdict), data = d[verdict == "no_FG_CS"], size = 2.4) +
    geom_text_repel(aes(label = pheno_short, colour = verdict),
                    size = 1.8, seed = 7, max.overlaps = 16, min.segment.length = 0.15,
                    segment.size = 0.2, segment.colour = "grey60", box.padding = 0.14, force = 1.1, max.time = 8) +
    facet_wrap(~ gene, ncol = 2) +
    scale_x_log10(limits = c(FLOOR * 0.7, 1.4), breaks = 10^seq(-4, 0), labels = c("1e-4", "1e-3", "0.01", "0.1", "1")) +
    scale_y_log10(limits = c(FLOOR * 0.7, 1.4), breaks = 10^seq(-4, 0), labels = c("1e-4", "1e-3", "0.01", "0.1", "1")) +
    scale_colour_manual(values = pal, name = NULL,
                        labels = c(confirmed = "confirmed (abf & SuSiE \u2265 0.8)", rescued = "SuSiE-only (\u2265 0.8)",
                                   reversed = "abf-only (overturned)", neither = "neither", no_FG_CS = "no FinnGen credible set"),
                        guide = guide_legend(override.aes = list(linetype = 0, size = 2.2, shape = c(19,19,19,19,0)))) +
    scale_shape_manual(values = c(confirmed = 19, rescued = 19, reversed = 19, neither = 19, no_FG_CS = 0), guide = "none") +
    scale_alpha_manual(values = c(`TRUE` = 0.35, `FALSE` = 1), guide = "none") +
    labs(x = "coloc.abf PP.H4 (single-signal)", y = "coloc.susie PP.H4 (max CS pair)",
         title = "SuSiE multi-signal pleiotropy review (chr12 candidates)",
         subtitle = "SuSiE confirms GDF11 pleiotropy but prunes SUOX hypothyroid/autoimmune; OLP coloc robust in both") +
    theme_journal(8) +
    theme(strip.text = element_text(face = "italic", size = 8.5),
          legend.position = "bottom", legend.text = element_text(size = 6.3),
          legend.key.size = unit(0.28, "cm"))
}

# Fig7B: pleiotropy heatmap
render_heatmap <- function() {
  PHENO_LAB <- c(K11_LICHEN_PLANUS_WIDE = "OLP (WIDE)", K11_ORAL_LICHEN_PLANUS_INCLAVO = "OLP (INCLAVO)",
                 K11_COELIAC = "Coeliac disease", T1D_WIDE = "Type 1 diabetes",
                 E4_HYTHY_AI_STRICT = "Hypothyroidism (autoimmune)", J10_NASALPOLYP = "Nasal polyps",
                 L12_LICHENSCLERATROPH = "Lichen sclerosus", E4_GRAVES_STRICT = "Graves disease",
                 RHEUMA_SEROPOS_STRICT = "RA (seropositive)", AUTOIMMUNE = "Any autoimmune disease")
  d <- fread("results/tables/pleiotropy_coloc_results.csv")
  cand <- c("GDF11", "SUOX", "EGFL8", "TAPBP", "HLA-DPA1")
  d <- d[gene %in% cand]
  d[, pheno_lab := factor(PHENO_LAB[pheno], levels = unname(PHENO_LAB))]
  d[, gene := factor(gene, levels = rev(cand))]
  d[, pph4 := pmin(PP.H4, 1)]
  ggplot(d, aes(x = pheno_lab, y = gene, fill = pph4)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", PP.H4), colour = ifelse(pph4 > 0.6, "white", "black")), size = 2.0) +
    scale_fill_gradientn(colours = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
                         values = c(0, 0.25, 0.5, 0.8, 1), limits = c(0, 1), name = "PP.H4") +
    scale_colour_identity() +
    labs(x = NULL, y = NULL,
         title = "Pleiotropy colocalisation matrix (coloc.abf PP.H4)",
         subtitle = "chr12 (GDF11/SUOX) pleiotropy is a shared causal variant; MHC-locus pleiotropy is not; OLP = internal control") +
    theme_journal(8) +
    theme(axis.text.x = element_text(angle = 38, hjust = 1, size = 6.3),
          axis.text.y = element_text(size = 7.5, face = "italic"),
          legend.position = "right", legend.key.width = unit(0.25, "cm"), legend.key.height = unit(0.5, "cm"),
          legend.title = element_text(size = 7))
}

# ============================================================
# Fig1: 设计(A 上半占满) + 下半 Venn(B 左) + funnel(C 右)
# ============================================================
fig1A <- render_design() + labs(tag = "A")
fig1C <- render_funnel() + labs(tag = "C") +
  labs(title = "Druggable genes with cis-eQTL support",
       subtitle = "5,801 of 9,261 druggable genes (62.6%)\nhave an eQTLGen cis-eQTL instrument")
fig1B_raw <- render_venn()
fig1B_tagged <- cowplot::ggdraw(fig1B_raw) + labs(tag = "B") + theme(plot.tag = element_text(face = "bold", size = 11))
bottom_row <- (fig1B_tagged | fig1C) + plot_layout(widths = c(1, 1.15))
fig1 <- (fig1A / bottom_row) + plot_layout(heights = c(1.15, 1)) &
  theme(plot.tag = element_text(face = "bold", size = 11))
save_main("fig01_design_sources", fig1, 180, 175)

# ============================================================
# Fig2-6: 复制已定稿 journal 版（重命名 fig02-fig06）
# ============================================================
jmap <- c("fig2_mr_results" = "fig02_mr_screens",
          "fig3_forest_tiers" = "fig03_candidates_forest",
          "fig4_coloc_evidence" = "fig04_colocalisation",
          "fig5_phewas_steiger" = "fig05_pleiotropy_steiger",
          "fig6_docking" = "fig06_docking")
for (src in names(jmap)) {
  for (ext in c("png", "pdf")) {
    f <- file.path("figures/journal", paste0(src, ".", ext))
    if (file.exists(f)) file.copy(f, file.path(MDIR, paste0(jmap[src], ".", ext)), overwrite = TRUE)
  }
}
cat("copied Fig2-6 from journal\n")

# ============================================================
# Fig7: SuSiE pleiotropy(A) + heatmap(B) —— 上下两面板
# ============================================================
fig7A <- render_susie_pleiotropy() + labs(tag = "A")
fig7B <- render_heatmap() + labs(tag = "B")
fig7 <- (fig7A / fig7B) + plot_layout(heights = c(1, 1)) &
  theme(plot.tag = element_text(face = "bold", size = 11))
save_main("fig07_pleiotropy_verdict", fig7, 180, 190)

cat("ALL MAIN FIGURES DONE\n")
