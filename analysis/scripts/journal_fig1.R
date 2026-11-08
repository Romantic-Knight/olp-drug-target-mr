# journal_fig1.R —— 论文主图 Fig1：研究设计示意（5 步漏斗）
# 每步一个独立列块：居中 x=1..5，宽度/高度随步骤递减，块间留固定间隙；
# 标题在块上方、数字在块内、说明在块下方。cohort 标注在 MR screen 上方高区。
# 输出 figures/journal/fig1_study_design.{png,pdf,tiff,eps}
suppressMessages({library(ggplot2)})
source("analysis/scripts/journal_theme.R")

# ---- 各步真实样本数/阈值 ----
n_druggable <- 9261; n_intersect <- 5801; n_tested <- 5040
n_wide <- 9526; n_incl <- 2491
n_coloc <- 44; n_candidates <- 5

# ---- 5 步：居中 x=1..5，宽度递减、高度递减 ----
xs <- 1:5
# 宽度递减但保证最小可读，且相邻块(x±width/2)不重叠：间距1，半宽和<1
widths  <- c(0.92, 0.80, 0.70, 0.60, 0.52)
heights <- c(0.95, 0.80, 0.66, 0.54, 0.46)
yc <- 2.05
cols <- c("#005BBB", "#1B9E77", "#2B6CB0", "#D1495B", "#C8102E")

fw <- data.frame(
  x = xs,
  xmin = xs - widths / 2, xmax = xs + widths / 2,
  ymin = yc - heights / 2, ymax = yc + heights / 2,
  fill = cols)

step_lab <- c("Druggable\ngenome", "eQTLGen\ncis-eQTL\noverlap", "MR screen\n(2 OLP\ndefinitions)",
              "Colocalisation\n(abf + SuSiE)", "Tiered\ncandidates")
steps <- data.frame(
  x = xs, step = step_lab,
  n = c(n_druggable, n_intersect, n_tested, n_coloc, n_candidates),
  sub = c("9,261 druggable\ngenes (4 sources)", "5,801 genes\nw/ cis-eQTL",
          "5,040 genes\ntested (Bonf 9.9e-6)", "44 abf + 44\nSuSiE runs",
          "5 candidates\n(1xT1,2xT2,2xT3)"))
steps <- cbind(steps, fw[, c("ymin", "ymax")])

# cohort 标注：MR screen (x=3) 上方高区
cohorts <- data.frame(x = c(2.74, 3.52), y = 3.98,
                      cohort = c("WIDE 9,526 cases", "INCLAVO 2,491 cases"))

# 底部数据源注释：x=1,2,3,4 间隔
ann <- data.frame(x = c(1, 2.4, 3.6, 4.75), y = 0.34,
                  txt = c("FinnGen R13\n(GRCh38)", "eQTLGen\n31,684", "1000G LD\n(EUR, n=503)",
                          "pQTL evidence\nINTERVAL/Knight/\ndeCODE 35,652"))

p <- ggplot() +
  geom_rect(data = fw, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = fw$fill, alpha = 0.30, colour = "grey30", linewidth = 0.35) +
  # 相邻块连线（同中线）
  geom_segment(data = data.frame(x = xs[-5] + widths[-5]/2, xend = xs[-1] - widths[-1]/2, y = yc, yend = yc),
               aes(x = x, xend = xend, y = y, yend = yend),
               colour = "grey55", linetype = "dashed", linewidth = 0.28) +
  # 步骤标题（块上方）
  geom_text(data = steps, aes(x = x, y = ymax + 0.55, label = step),
            hjust = 0.5, size = 2.8, fontface = "bold", lineheight = 0.85) +
  # 块内大数字
  geom_text(data = steps, aes(x = x, y = ymax - 0.05, label = n),
            hjust = 0.5, vjust = 1, size = 3.2, fontface = "bold", colour = "black") +
  # 块下方小字说明
  geom_text(data = steps, aes(x = x, y = ymin - 0.05, label = sub),
            hjust = 0.5, vjust = 1, size = 2.1, colour = "grey25", lineheight = 0.82) +
  # cohort 标注（虚线引到 MR screen 块顶）
  geom_segment(data = cohorts, aes(x = x, xend = x, y = 3.98 - 0.13, yend = 2.55),
               linetype = "dashed", colour = "grey55", linewidth = 0.26) +
  geom_text(data = cohorts, aes(x = x, y = y, label = cohort),
            hjust = 0.5, size = 2.4, fontface = "bold", colour = "grey15") +
  # 底部注释
  geom_text(data = ann, aes(x = x, y = y, label = txt),
            hjust = 0.5, vjust = 1, size = 2.1, colour = "grey30", lineheight = 0.8) +
  coord_cartesian(xlim = c(0.28, 5.72), ylim = c(-0.02, 4.32)) +
  theme_journal(8) +
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank(),
        plot.margin = margin(6, 6, 6, 6)) +
  labs(title = "Study design: druggable-genome MR of oral lichen planus",
       subtitle = "Cis-eQTL instruments (eQTLGen) -> FinnGen R13 OLP (two case definitions) -> coloc.abf/SuSiE -> pQTL & docking")

save_journal("fig1_study_design", p, 180, 105)
