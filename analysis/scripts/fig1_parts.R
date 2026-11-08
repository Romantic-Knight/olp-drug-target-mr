# fig1_parts.R —— 主图 Fig1 的三个面板渲染器（单一来源）
#   render_design()  设计漏斗（5 步）
#   render_venn()    四来源可成药基因 Venn
#   render_funnel()  可成药 ∩ eQTLGen 支持柱
# 由 build_final_fig01.R 与 build_main_figures.R 共同 source，避免两份副本漂移。
# 需在项目根目录下运行（fread 使用项目相对路径）。
if (!exists("theme_journal")) source("analysis/scripts/journal_theme.R")
suppressMessages({library(ggplot2); library(data.table); library(VennDiagram); library(cowplot)})

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

# Fig1B: druggable Venn
# B 与 C 并排时必须"同构"：两者都是 ggplot，都自带 title/subtitle 和一块 panel。
# 若 B 用 cowplot::ggdraw() 且不设标题，它的 panel 会占满整幅图，与 C（标题在上、
# 坐标轴在下，panel 较小）并排时组合器对齐的是两者的 panel，B 于是被压到 C 的
# panel 位置而整体下沉（实测标题比 C 低约 120 px）。这里给 ggdraw 补上同样的
# title/subtitle（并恢复被 theme_nothing 抹掉的 plot.title/subtitle 元素），
# 使两者的 panel 处于同一水平线；draw_grob 采用 0,0,1,1 让 Venn 铺满该 panel。
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
    cowplot::draw_grob(grob, x = 0, y = 0, width = 1, height = 1) +
    ggplot2::labs(title = "Druggable genes: four-source union",
                  subtitle = sprintf("Pharos / DGIdb / Finan 2017 / OpenTargets (union n = %s)",
                                     format(nrow(fin), big.mark = ","))) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 9, hjust = 0.5,
                                         colour = "black", margin = ggplot2::margin(b = 2)),
      plot.subtitle = ggplot2::element_text(size = 7, hjust = 0.5,
                                            colour = "grey35", margin = ggplot2::margin(b = 2)),
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    )
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
