# build_final_mr_figures.R —— 最终 8 主图方案中的 MR 5 张（Fig1-5）
# 映射：
#   Fig1 = 原 Fig1（设计 + Venn + eQTL 柱）      -> 保留 main_figures/fig01_design_sources
#   Fig2 = 原 Fig2（WIDE/INCLAVO Manhattan+volcano）-> 保留 fig02_mr_screens
#   Fig3 = 原 Fig3（候选森林图）+ 原 Fig6（分子对接） -> 合并 fig03_candidates_docking
#   Fig4 = 原 Fig4（共定位）                        -> 保留 fig04_colocalisation
#   Fig5 = 原 Fig5（PheWAS + Steiger）+ 原 Fig7（SuSiE + PP.H4 矩阵）-> 合并 fig05_pleiotropy_susie
# Fig1/2/4 已在 main_figures 中定稿，本脚本只重渲染 Fig3（合并）与 Fig5（合并）。
suppressMessages({library(ggplot2); library(data.table); library(ggrepel); library(VennDiagram); library(grid); library(patchwork); library(cowplot); library(scales)})
source("analysis/scripts/journal_theme.R")

MDIR <- "figures/main_figures"; dir.create(MDIR, recursive = TRUE, showWarnings = FALSE)
save_main <- function(base, plot, width_mm, height_mm) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  ggplot2::ggsave(file.path(MDIR, paste0(base, ".png")), plot, width = w, height = h, dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(MDIR, paste0(base, ".pdf")), plot, width = w, height = h, device = grDevices::cairo_pdf)
  cat("saved", base, "\n")
}
TAG <- ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 11))

# 在独立环境中载入单图 journal 脚本并捕获绘图对象（save_journal 会重写 figures/journal，无害）
load_plot <- function(script, name) {
  e <- new.env(parent = globalenv())
  sys.source(script, envir = e)
  e[[name]]
}

# ============================================================
# Fig3 = 原 Fig3（森林图 A）+ 原 Fig6（分子对接 B）
# ============================================================
forest  <- load_plot("analysis/scripts/journal_fig3.R", "p")
docking <- load_plot("analysis/scripts/journal_fig6.R", "p")
# 修复 GDF11 条带标签在右缘被裁剪的问题：右侧条带标签改为水平并留出右缘余量
docking <- docking + ggplot2::theme(
  strip.text.y = ggplot2::element_text(angle = 0, hjust = 0.5),
  plot.margin = ggplot2::margin(6, 16, 6, 6))
fig3 <- (forest + labs(tag = "A")) / (docking + labs(tag = "B")) +
  patchwork::plot_layout(heights = c(1, 1.05)) & TAG
save_main("fig03_candidates_docking", fig3, 180, 215)

# ============================================================
# Fig5 = 原 Fig5（PheWAS A + Steiger B）+ 原 Fig7（SuSiE C + 矩阵 D）
# ============================================================
e5 <- new.env(parent = globalenv())
sys.source("analysis/scripts/journal_fig5.R", envir = e5)
phewas  <- e5$pa
steiger <- e5$pb

# --- 原 Fig7A：SuSiE 多信号裁决（复刻 build_main_figures.R）---
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
                        guide = guide_legend(nrow = 2, byrow = TRUE,
                                               override.aes = list(linetype = 0, size = 2.2, shape = c(19,19,19,19,0)))) +
    scale_shape_manual(values = c(confirmed = 19, rescued = 19, reversed = 19, neither = 19, no_FG_CS = 0), guide = "none") +
    scale_alpha_manual(values = c(`TRUE` = 0.35, `FALSE` = 1), guide = "none") +
    labs(x = "coloc.abf PP.H4 (single-signal)", y = "coloc.susie PP.H4 (max CS pair)",
         title = "SuSiE multi-signal pleiotropy review\n(chr12 candidates)",
         subtitle = "SuSiE confirms GDF11 pleiotropy but prunes SUOX\nhypothyroid/autoimmune; OLP coloc robust in both") +
    theme_journal(8) +
    theme(strip.text = element_text(face = "italic", size = 8.5),
          legend.position = "bottom", legend.text = element_text(size = 6.0),
          legend.key.size = unit(0.28, "cm"),
          legend.box = "horizontal",
          legend.margin = margin(t = -2, b = 0))
}

# --- 原 Fig7B：多效性共定位矩阵（复刻 build_main_figures.R）---
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
         title = "Pleiotropy colocalisation matrix\n(coloc.abf PP.H4)",
         subtitle = "chr12 (GDF11/SUOX) pleiotropy is a shared causal variant;\nMHC-locus pleiotropy is not; OLP = internal control") +
    theme_journal(8) +
    theme(axis.text.x = element_text(angle = 38, hjust = 1, size = 6.3),
          axis.text.y = element_text(size = 7.5, face = "italic"),
          legend.position = "right", legend.key.width = unit(0.25, "cm"), legend.key.height = unit(0.5, "cm"),
          legend.title = element_text(size = 7))
}

susie   <- render_susie_pleiotropy()
heatmap <- render_heatmap()

# 2x2 拼接：上排 A PheWAS | B Steiger；下排 C SuSiE | D 矩阵
# 画布加宽（220 mm，约定双栏 180 mm 之上的宽版），给分面内的表型标签留出更多
# 横向空间，避免挤在一起。
# 面板原为整幅 180 mm 宽设计，标题/副标题/图例按整幅宽度排布；压到半幅
# （约 86 mm）后会出现标题跨列重叠与图例被裁切，故在此把标题与副标题手工
# 换行，并让 C 的底图例分两行。
phewas <- phewas + labs(
  title = "PheWAS: all 5 candidate IVs are pleiotropic\nacross FinnGen R13 phenotypes",
  subtitle = paste0("dashed = genome-wide 5e-8; red diamond = OLP endpoint;\n",
                    "same/opposite relative to each gene's own OLP effect;\n",
                    "usable phenotypes per IV: 415/308/392/298/310"))
# (B 的标题/副标题在下面重建 steiger 时一并设置：此处若赋值会被后续
#  steiger <- ggplot(...) 覆盖，故 B 面板不能在这里定义标题。)

# A 在 2x2 中只有约 86 mm 宽（5 个分面各约 15 mm），表型标注会被分面边界裁掉。
# 关闭 panel 裁剪让标签可以伸到分面外，并加宽分面间距腾出空间；同时略缩小标注、
# 收紧 max.overlaps 以减少分面内叠字。只改最终图对象，不动 journal_fig5.R。
phewas <- phewas +
  coord_cartesian(clip = "off") +
  theme(panel.spacing.x = unit(7, "pt"))
for (i in seq_along(phewas$layers)) {
  if (inherits(phewas$layers[[i]]$geom, "GeomTextRepel")) {
    phewas$layers[[i]]$aes_params$size <- 1.35
    phewas$layers[[i]]$geom_params$max.overlaps <- 9
    # 分面内密集标签需较大排斥力 + 足量迭代才会收敛到不重叠
    phewas$layers[[i]]$geom_params$force <- 2
    phewas$layers[[i]]$geom_params$force_pull <- 0.5
    phewas$layers[[i]]$geom_params$max.time <- 10
    phewas$layers[[i]]$geom_params$max.iter <- 100000
    phewas$layers[[i]]$geom_params$box.padding <- 0.1
    phewas$layers[[i]]$geom_params$point.padding <- 0.1
    phewas$layers[[i]]$geom_params$seed <- 1   # 固定随机起点，保证可复现
  }
}

# B：同一基因在 WIDE 与 INCLAVO 两个方向上各有一点，会画出 8 点，与
# "all 5 candidates" 不符。只保留 WIDE 方向（标题/轴标签本就写 WIDE），
# 每基因恰一点。
stc <- fread("results/tables/steiger_direction_test.csv")
stw <- stc[grepl("LICHEN_PLANUS", pheno)]
cand <- c("GDF11","SUOX","EGFL8","TAPBP","HLA-DPA1")
stw[, is_cand := ifelse(gene %in% cand, gene, "")]
stw <- stw[pheno == "K11_LICHEN_PLANUS_WIDE" | is_cand == ""]
steiger <- ggplot(stw, aes(x = z_exp, y = z_out)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(aes(colour = is_cand == ""), size = 1.3, alpha = 0.8) +
  geom_point(data = stw[is_cand != ""], aes(x = z_exp, y = z_out),
             shape = 18, size = 3, colour = "#C8102E") +
  ggrepel::geom_text_repel(data = stw[is_cand != ""], aes(label = is_cand),
             size = 2.6, seed = 1, colour = "#C8102E") +
  scale_colour_manual(values = c(`TRUE` = "#9AA5B1", `FALSE` = "#C8102E"), guide = "none") +
  labs(x = "eQTLGen cis-eQTL Z (expression)", y = "FinnGen OLP WIDE MR Z (outcome)",
       title = "Steiger direction test:\nforward orientation for 4 of 5 candidates",
       subtitle = "red = tiered candidate; EGFL8 forward but P = 0.082;\nGDF11 does not pass in INCLAVO (P = 0.72)") +
  theme_journal(8)

fig5 <- ((phewas + labs(tag = "A")) | (steiger + labs(tag = "B"))) /
        ((susie + labs(tag = "C")) | (heatmap + labs(tag = "D")))
fig5 <- fig5 + patchwork::plot_layout(heights = c(1, 1)) & TAG
save_main("fig05_pleiotropy_susie", fig5, 220, 220)

cat("MR FINAL FIGURES DONE\n")
