# 分级靶点 MR 森林图 —— 口腔扁平苔藓（OLP）drug-target MR Step 1-6 汇总
# 输入: results/tables/mr_olp_hits_crosscohort.csv（22 个 WIDE Bonf hits 双队列 MR）
#       results/tables/coloc_olp_results.csv（abf PP.H4 注释）
# 输出: figures/fig_mr_forest_tiers.png / .pdf
suppressMessages({library(ggplot2)})

hits <- read.csv("results/tables/mr_olp_hits_crosscohort.csv", check.names = FALSE)
coloc <- read.csv("results/tables/coloc_olp_results.csv", check.names = FALSE)

# 最终定级（Step 5c 定稿；SUOX 因 chr12 LD 非独立降级 TIER2，EGFL8 维持 TIER2）
final_tier <- c(GDF11 = "TIER1", SUOX = "TIER2", EGFL8 = "TIER2",
                TAPBP = "TIER3", `HLA-DPA1` = "TIER3")
genes <- names(final_tier)

# 从 22 hits 表中取 5 个定级靶点的行（应全部存在）
sub <- hits[hits$gene %in% genes, c("gene","Z_MR_wide","p_wide","beta_MR_wide","OR_wide",
                                    "Z_MR_incl","p_incl","beta_MR_incl","OR_incl")]
stopifnot(nrow(sub) == length(genes))

mk <- function(g, cohort, or, beta, z, p) {
  se <- abs(beta / z)                       # 单 IV Wald: se = |beta|/|Z|
  data.frame(gene = g, cohort = cohort, or = or,
             lo = exp(beta - 1.96 * se), hi = exp(beta + 1.96 * se),
             p = p, stringsAsFactors = FALSE)
}
long <- do.call(rbind, lapply(seq_len(nrow(sub)), function(i) {
  r <- sub[i, ]
  rbind(mk(r$gene, "WIDE",   r$OR_wide,   r$beta_MR_wide,   r$Z_MR_wide,   r$p_wide),
        mk(r$gene, "INCLAVO", r$OR_incl,  r$beta_MR_incl,  r$Z_MR_incl,   r$p_incl))
}))
long$tier <- factor(final_tier[long$gene], levels = c("TIER1","TIER2","TIER3"))
long$cohort <- factor(long$cohort, levels = c("WIDE","INCLAVO"))
# y 轴顺序：TIER1 在顶部；每基因 WIDE 在 INCLAVO 上方
ord <- unlist(lapply(genes, function(g) paste(g, c("WIDE","INCLAVO"))))
long$lbl <- factor(paste(long$gene, long$cohort), levels = rev(ord))

# abf coloc PP.H4 注释
pp <- coloc[coloc$gene %in% genes & coloc$cohort %in% c("WIDE","INCLAVO"),
            c("gene","cohort","PP.H4")]
pp$lbl <- paste(pp$gene, pp$cohort)
pp$txt <- sprintf("PP.H4=%.2f", pp$PP.H4)
pp <- pp[match(levels(long$lbl), pp$lbl), ]
pp$tier <- final_tier[pp$gene]
pp$cohort <- factor(pp$cohort, levels = c("WIDE","INCLAVO"))

long$ptxt <- ifelse(long$p < 1e-6, sprintf("P=%.1e", long$p),
                    ifelse(long$p < 0.001, sprintf("P=%.2e", long$p),
                           sprintf("P=%.4g", long$p)))
# 数值文本放于最右注释列的位置 x（超过全部 CI 上界）
ann_x <- max(long$hi) * 1.35
xmax  <- max(long$hi) * 1.75

tier_col <- c(TIER1 = "#C8102E", TIER2 = "#005BBB", TIER3 = "#7F7F7F")
shape_col <- c(WIDE = 19, INCLAVO = 1)

p <- ggplot(long, aes(x = or, y = lbl)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = tier), height = 0.18, linewidth = 0.7) +
  geom_point(aes(shape = cohort, colour = tier), size = 2.6) +
  geom_text(aes(x = or, label = ptxt), hjust = -0.35, vjust = -0.9, size = 2.6, colour = "grey25") +
  geom_text(data = pp, aes(x = ann_x, y = lbl, label = txt),
            hjust = 0, size = 2.6, colour = "grey30", fontface = "italic") +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4, 8, 16, 32), limits = c(0.1, xmax)) +
  scale_colour_manual(values = tier_col, name = "Final tier") +
  scale_shape_manual(values = shape_col, name = "OLP cohort") +
  labs(x = "MR Wald-ratio OR (eQTLGen cis-eQTL \u2192 FinnGen R13 OLP; log scale)",
       y = NULL,
       title = "Drug-target MR estimates of tiered candidate targets \u2014 oral lichen planus",
       subtitle = "Single top cis-eQTL IV per gene; error bars = 95% CI; dashed line OR = 1; italic right column = coloc.abf PP.H4") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey35", size = 9))
ggsave("figures/fig_mr_forest_tiers.png", p, width = 9.2, height = 5.6, dpi = 300)
ggsave("figures/fig_mr_forest_tiers.pdf", p, width = 9.2, height = 5.6, device = cairo_pdf)
cat("rows:", nrow(long), "| genes:", paste(genes, collapse = ","), "\n")
cat("GDF11 INCLAVO OR/CI:", subset(long, gene == "GDF11" & cohort == "INCLAVO")$or,
    "| lo/hi:", subset(long, gene == "GDF11" & cohort == "INCLAVO")$lo,
    subset(long, gene == "GDF11" & cohort == "INCLAVO")$hi, "\n")
