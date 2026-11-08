# journal_fig4.R —— 论文主图 Fig4：共定位证据
# a) coloc.abf PP.H4（22 基因×2 队列，PP.H4>0.8 strong 标色）
# b) coloc.susie vs coloc.abf PP.H4（多信号复核）
suppressMessages({library(ggplot2); library(data.table)})
source("analysis/scripts/journal_theme.R")

abf <- fread("results/tables/coloc_olp_results.csv")
sus <- fread("results/tables/susie_coloc_results.csv")[, .(gene, cohort, pp_susie = pp_h4_max)]
mr  <- fread("results/tables/mr_olp_wide_results.csv")[significant_bonf == TRUE, .(gene, p_wide = p_MR)]
d   <- merge(abf, mr, by = "gene", all.x = TRUE)
d[, strong := PP.H4 > 0.8]
# 22 WIDE hits 排序（按 PP.H4 升序，便于条图）
hitgenes <- unique(abf[cohort == "WIDE"][order(-PP.H4)]$gene)
d[, gene := factor(gene, levels = rev(hitgenes))]
d[, cohort := factor(cohort, levels = c("WIDE", "INCLAVO"))]

# a) abf PP.H4 条形（双队列并排）
pa <- ggplot(d, aes(x = gene, y = PP.H4, fill = strong)) +
  geom_col(width = 0.72, alpha = 0.92) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey40", linewidth = 0.35) +
  facet_grid(. ~ cohort, scales = "free_x") +
  scale_fill_manual(values = c(`FALSE` = "#9AA5B1", `TRUE` = "#D1495B"),
                    name = "PP.H4 > 0.8", labels = c("no", "yes")) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = "coloc.abf PP.H4 (shared causal variant)",
       title = "Colocalisation (abf): 22 WIDE Bonferroni hits",
       subtitle = "dashed line = PP.H4 0.8; WIDE ordered by PP.H4") +
  theme_journal(8) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 5.6))

# b) susie vs abf PP.H4（双队列点，PP.H4>0.8 参考线）
db <- merge(sus, abf[, .(gene, cohort, pp_abf = PP.H4)], by = c("gene", "cohort"))
db[, tier_keep := fifelse(pp_susie >= 0.8 & pp_abf >= 0.8, "both",
                    fifelse(pp_susie >= 0.8, "susie", "neither"))]
pb <- ggplot(db, aes(x = pp_abf, y = pp_susie, colour = tier_keep)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey65", linewidth = 0.3) +
  geom_vline(xintercept = 0.8, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  geom_point(size = 1.8, alpha = 0.85) +
  facet_wrap(~ cohort, nrow = 1) +
  scale_colour_manual(values = c(both = "#D1495B", susie = "#2B6CB0", neither = "#9AA5B1"),
                      name = NULL,
                      labels = c(both = "abf + SuSiE \u2265 0.8", susie = "SuSiE only \u2265 0.8",
                                 neither = "neither \u2265 0.8")) +
  scale_x_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25)) +
  labs(x = "coloc.abf PP.H4", y = "coloc.susie PP.H4 (max CS pair)",
       title = "Multi-signal (SuSiE) vs single-signal (abf) colocalisation",
       subtitle = "dotted line = 1:1 diagonal (points above gain from multi-CS); dashed = PP.H4 0.8; 6/44 runs with no credible set omitted") +
  theme_journal(8)

suppressMessages(library(patchwork))
fig4 <- (pa + labs(tag = "A")) / (pb + labs(tag = "B")) &
  ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 11))
save_journal("fig4_coloc_evidence", fig4, 180, 150)
