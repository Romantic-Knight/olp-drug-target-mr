# journal_fig3.R —— 论文主图 Fig3：候选靶点 MR 森林图（TIER1-3，双队列）
# 数据：final_target_evidence.csv（已含 OR/CI/P/susie PP.H4/verdict）——单一权威来源
suppressMessages({library(ggplot2); library(data.table)})
source("analysis/scripts/journal_theme.R")

d <- fread("results/tables/final_target_evidence.csv")
tier_lev <- c("TIER1", "TIER2", "TIER3")
d[, tier := factor(tier, levels = tier_lev)]

# 长表：每基因 2 行（WIDE / INCLAVO），CI 字符串拆分
parse_ci <- function(s) as.numeric(strsplit(s, "-", fixed = TRUE)[[1]])
long <- rbind(
  data.table(gene = d$gene, tier = d$tier, cohort = "WIDE",
             or = d$or_wide, lo = sapply(d$ci_wide, function(s) parse_ci(s)[1]),
             hi = sapply(d$ci_wide, function(s) parse_ci(s)[2]),
             p = d$p_wide, susie = d$susie_pph4_wide),
  data.table(gene = d$gene, tier = d$tier, cohort = "INCLAVO",
             or = d$or_incl, lo = sapply(d$ci_incl, function(s) parse_ci(s)[1]),
             hi = sapply(d$ci_incl, function(s) parse_ci(s)[2]),
             p = d$p_incl, susie = d$susie_pph4_incl))
long[, cohort := factor(cohort, levels = c("WIDE", "INCLAVO"))]
# y 轴：TIER1 在上，每基因 WIDE 在上 INCLAVO 在下
gorder <- c("GDF11", "SUOX", "EGFL8", "TAPBP", "HLA-DPA1")
ord <- unlist(lapply(rev(gorder), function(g) paste(g, c("INCLAVO", "WIDE"))))
long[, lbl := factor(paste(gene, cohort), levels = ord)]
long[, ptxt := fifelse(p < 1e-6, sprintf("P = %.1e", p),
                 fifelse(p < 0.001, sprintf("P = %.2e", p), sprintf("P = %.4g", p)))]

p <- ggplot(long, aes(x = or, y = lbl)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = tier), height = 0.16, linewidth = 0.6) +
  geom_point(aes(shape = cohort, colour = tier), size = 2.4) +
  geom_text(aes(x = or, label = ptxt), hjust = -0.3, vjust = -1.2, size = 2.2, colour = "grey25") +
  # PP.H4 注释：右对齐放在面板右缘（x 取对数坐标上限内侧）
  geom_text(aes(x = max(long$hi) * 1.9 * 0.97, y = lbl, label = sprintf("PP.H4 %.2f", susie)),
            hjust = 1, size = 2.0, colour = "grey30", fontface = "italic") +
  scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8, 16, 32)) +
  scale_colour_manual(values = tier_colours, name = "Final tier") +
  scale_shape_manual(values = c(WIDE = 19, INCLAVO = 1), name = "OLP cohort") +
  # the left limit was hard-coded at 0.12; take the range from the data so that
  # estimates with a lower bound below it are not clipped
  coord_cartesian(xlim = c(min(long$lo) * 0.85, max(long$hi) * 1.9), clip = "off") +
  labs(x = "MR Wald-ratio OR (eQTLGen cis-eQTL \u2192 FinnGen R13 OLP; log scale)",
       y = NULL,
       title = "Tiered candidate drug targets for oral lichen planus",
       subtitle = paste0("Single top cis-eQTL IV per gene; error bars = 95% CI; dashed line OR = 1; ",
                         "OR < 1 = higher expression associated with lower risk; ",
                         "right column = coloc.susie PP.H4")) +
  theme_journal(8) +
  theme(plot.margin = margin(6, 30, 6, 6))

save_journal("fig3_forest_tiers", p, 180, 100)
