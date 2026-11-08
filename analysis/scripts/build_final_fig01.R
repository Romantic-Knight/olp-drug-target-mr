# build_final_fig01.R —— 最终主图 Fig1（只写 fig01_design_sources，不动其它主图）
# 版面：上 A 设计漏斗（整宽）；下排 B 可成药 Venn | C eQTL 支持柱。
#
# B/C 对齐修复：原实现里 B 是 cowplot::ggdraw（panel 占满整幅图）、C 是普通
# ggplot（panel 之上有标题、之下有 x 轴）。并排时组合器对齐的是两者的 **panel**
# 区域，B 因此被压到 C 的 panel 位置而整体下沉，标题比 C 低约 120 px，
# 看起来 B 与 C 不在同一水平线上。修法是把 B 也做成同构的 ggplot（Venn grob
# 用 annotation_custom 铺满 panel、标题交给 labs()），再用 patchwork 并排，
# 两者标题/副标题即天然对齐。
suppressMessages({library(ggplot2); library(patchwork)})
source("analysis/scripts/journal_theme.R")
source("analysis/scripts/fig1_parts.R")

MDIR <- "figures/main_figures"
TAG <- ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 11))

fig1A <- render_design() + labs(tag = "A")
fig1C <- render_funnel() + labs(tag = "C") +
  labs(title = "Druggable genes with cis-eQTL support",
       subtitle = "5,801 of 9,261 druggable genes (62.6%)\nhave an eQTLGen cis-eQTL instrument")
fig1B <- render_venn() + labs(tag = "B")

# A 整宽在上；下排 B | C。三个都是同构的 ggplot，patchwork 对齐 panel 后
# B 与 C 的标题/副标题天然落在同一水平线。
bottom_row <- (fig1B | fig1C) + patchwork::plot_layout(widths = c(1, 1.15))
fig1 <- (fig1A / bottom_row) + patchwork::plot_layout(heights = c(1.15, 1)) & TAG

w <- 180 / 25.4; h <- 175 / 25.4
ggplot2::ggsave(file.path(MDIR, "fig01_design_sources.png"), fig1,
                width = w, height = h, dpi = 300, bg = "white")
ggplot2::ggsave(file.path(MDIR, "fig01_design_sources.pdf"), fig1,
                width = w, height = h, device = grDevices::cairo_pdf)
cat("saved fig01_design_sources\n")
