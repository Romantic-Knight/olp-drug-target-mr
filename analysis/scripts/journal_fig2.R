# journal_fig2.R —— 论文主图 Fig2：双队列 MR 结果
# a) WIDE Manhattan  b) INCLAVO Manhattan  c) WIDE Volcano
# 数据来源 mr_olp_wide/inclavo_results.csv（含 chr/pos/Z_MR/p_MR/OR/显著标记）
suppressMessages({library(ggplot2); library(data.table); library(ggrepel)})
source("analysis/scripts/journal_theme.R")

# 字号整体下调一档（base 8 -> 7.5；标题/刻度/基因标注同比缩小），并把画布
# 由 180x140 mm 增至 180x148 mm。宽度仍是期刊双栏上限 180 mm，只加高，
# 使文字相对图幅变小并给 chr6 附近的密集基因标注留出排斥空间。
FS <- 7.5

w  <- fread("results/tables/mr_olp_wide_results.csv")
i  <- fread("results/tables/mr_olp_inclavo_results.csv")
cc <- fread("results/tables/mr_olp_hits_crosscohort.csv")

# ---- 染色体位置累积（GRCh37 连续坐标）----
prep <- function(d) {
  tmp <- d[, .(maxpos = max(pos)), by = chr][order(chr)]
  tmp[, maxpos := as.numeric(maxpos)]
  tmp[, off := cumsum(shift(maxpos, fill = 0))]
  tmp[, off := as.numeric(off)][, mid := off + maxpos / 2]
  d <- merge(d, tmp[, .(chr, off, mid)], by = "chr")
  d[, xpos := pos + off]
  d
}
w <- prep(w); i <- prep(i)

# 显著标记：WIDE Bonf；双队列 Bonf 复制的基因（跨队列复制）
bonf_w <- 0.05 / 5040
w[, repl_class := fifelse(significant_bonf & gene %in% cc[repl_bonf == TRUE]$gene, "repl",
                    fifelse(significant_bonf, "bonf", "ns"))]
i[, repl_class := fifelse(significant_bonf, "bonf", "ns")]
w[, label := fifelse(significant_bonf, gene, "")]
i[, label := fifelse(significant_bonf, gene, "")]

manh <- function(d, repl = TRUE, title, subt) {
  pal <- if (repl) c(ns = "#B0B0B0", bonf = "#2B6CB0", repl = "#D1495B") else c(ns = "#B0B0B0", bonf = "#2B6CB0")
  # 染色体刻度间隔显示（1,3,5,...,21），避免 18-22 密集重叠；21/22 mid 过近，仅标 21
  keep <- d$chr %% 2 == 1
  p <- ggplot(d, aes(x = xpos, y = -log10(p_MR))) +
    geom_point(aes(colour = repl_class), size = 0.9, alpha = 0.8) +
    geom_hline(yintercept = -log10(bonf_w), linetype = "dashed", colour = "grey40", linewidth = 0.35) +
    geom_text_repel(aes(label = label, colour = repl_class), size = 1.9, max.overlaps = 60,
                    seed = 1, min.segment.length = 0.3) +
    scale_colour_manual(values = pal, guide = "none") +
    scale_x_continuous(breaks = d$mid[keep], labels = d$chr[keep], expand = expansion(mult = c(0.005, 0.01))) +
    labs(x = "Chromosome (top cis-eQTL position, GRCh37)", y = expression(-log[10](italic(P))),
         title = title, subtitle = subt) +
    theme_journal(FS) +
    theme(axis.text.x = element_text(size = 6.1), panel.grid.major = element_blank(),
          plot.title = element_text(size = 7.8), plot.margin = margin(6, 4, 6, 6))
  p
}

# volcano：x=Z_MR, y=-log10(p_MR)，语义零点 x=0
w[, lab_v := fifelse(significant_bonf, gene, "")]
vol <- ggplot(w, aes(x = Z_MR, y = -log10(p_MR))) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = -log10(bonf_w), linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_point(aes(colour = repl_class), alpha = 0.75) +
  geom_text_repel(aes(label = lab_v, colour = repl_class), size = 1.9, max.overlaps = 50,
                  seed = 1, min.segment.length = 0.3) +
  scale_colour_manual(values = c(ns = "#B0B0B0", bonf = "#2B6CB0", repl = "#D1495B"),
                      name = "WIDE Bonferroni hit",
                      labels = c(ns = "not significant", bonf = "WIDE Bonf.", repl = "replicated in\nINCLAVO (Bonf.)")) +
  scale_x_continuous(limits = c(min(w$Z_MR) - 0.5, max(w$Z_MR) + 0.5)) +
  labs(x = "MR Wald-ratio Z (per-SD expression effect on OLP)",
       y = expression(-log[10](italic(P))),
       title = "Volcano: WIDE MR estimates",
       subtitle = paste0("n = 5,040 genes; WIDE Bonferroni hits highlighted")) +
  theme_journal(FS) +
  theme(legend.position = "inside", legend.position.inside = c(0.02, 0.5),
        legend.justification = c(0, 0.5),
        legend.background = element_rect(colour = "grey75", fill = alpha("white", 0.85), linewidth = 0.3),
        legend.margin = margin(2, 3, 2, 3),
        plot.margin = margin(6, 10, 6, 6))

# 面板字母统一：11 pt 粗体，与 Fig1/3/5 及本图 C 一致。
# 修复记录：原写法把 panel tag 主题放在组合式末尾
#   `(png_a / png_b) | vol + plot_layout(...) & theme(...)`
# 而 R 的运算符优先级是 `+` > `&` > `|`，因此该 `&` 只作用到
# `vol + plot_layout(...)` 这一支；嵌套的 `(png_a / png_b)` 拿不到它，
# 退回 ggplot2 默认的 `plot.tag = rel(1.2)`（base 7.5 -> 9 pt 且非粗体）。
# 结果就是 C 明显比 A、B 粗且大（实测 cap 高 32 px vs 27 px @300 dpi）。
# 现在把 tag 主题显式加到三个组件图各自身上，不再依赖 `&` 的传播范围。
TAG <- ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 11))
png_a <- manh(w, TRUE,  "Manhattan: WIDE (n = 9,526)",
              paste0("n = 5,040 genes; dashed line = Bonferroni (P = ", format(bonf_w, scientific = TRUE, digits = 2), ")")) +
  labs(tag = "A") + TAG
png_b <- manh(i, FALSE, "Manhattan: INCLAVO (n = 2,491)",
              paste0("n = 5,040 genes; dashed line = Bonferroni (P = ", format(bonf_w, scientific = TRUE, digits = 2), ")")) +
  labs(tag = "B") + TAG
vol <- vol + labs(tag = "C") + TAG

# 布局：a/b 上下并排（左）、c 在右侧（等高）。宽度保持本图此前实际渲染出的
# 1:1——历史版本里 `plot_layout(widths = c(1.9, 1))` 被吞掉了（见上），
# 实际成稿一直是 1:1，本轮只修 tag，不动版面比例。
suppressMessages(library(patchwork))
fig2 <- ((png_a / png_b) | vol) +
  patchwork::plot_layout(widths = c(1, 1))
save_journal("fig2_mr_results", fig2, 180, 148)
