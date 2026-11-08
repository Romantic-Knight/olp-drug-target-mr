# journal_fig6.R —— 论文主图 Fig6：分子对接（SUOX 8 good + GDF11）
# 数据：results/tables/docking_results.csv（target/ligand/best_affinity_kcal_mol）
suppressMessages({library(ggplot2); library(data.table)})
source("analysis/scripts/journal_theme.R")

d <- fread("results/tables/docking_results.csv")
d[, good := best_affinity_kcal_mol < -7]
# 只显示 good（< -7）命中；排序按亲和力
d <- d[good == TRUE][order(best_affinity_kcal_mol)]
d[, ligand := factor(ligand, levels = rev(ligand))]
d[, target_lab := fifelse(target == "SUOX", "SUOX (2A9A, apo)", "GDF11 (5E4G)")]

p <- ggplot(d, aes(x = ligand, y = best_affinity_kcal_mol, fill = target)) +
  geom_col(width = 0.66) +
  coord_flip() +
  geom_hline(yintercept = -7, linetype = "dashed", colour = "grey35", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", best_affinity_kcal_mol)),
            hjust = -0.15, size = 2.4, colour = "grey20") +
  facet_grid(target_lab ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c(SUOX = "#2B6CB0", GDF11 = "#D1495B"), guide = "none") +
  scale_y_continuous(limits = c(min(d$best_affinity_kcal_mol) - 0.8, 0.5),
                     breaks = seq(-9, -6.5, 0.5)) +
  labs(x = NULL, y = "AutoDock Vina affinity (kcal/mol)",
       title = "Molecular docking of druggable hits (threshold < -7 kcal/mol)",
       subtitle = "dashed line = good-docking threshold; SUOX receptor 2A9A chain A (apo); GDF11 5E4G") +
  theme_journal(8)

save_journal("fig6_docking", p, 180, 110)
