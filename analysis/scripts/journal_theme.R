# journal_theme.R —— 期刊统一主题模块（figure-style 规范）
# 字号阶梯 8/7/6（base=8 用于标题/轴标签，7 图例/注释，6 ticks）
# 外刻 ticks、无框 legend、色盲安全（同前图），Arial(Helvetica) 优先
# 用法: source("analysis/scripts/journal_theme.R")

theme_journal <- function(base = 8) {
  ggplot2::theme_classic(base_size = base, base_family = "sans") +
    ggplot2::theme(
      axis.ticks.length = ggplot2::unit(1.6, "pt"),
      axis.line = ggplot2::element_line(colour = "grey25", linewidth = 0.4),
      axis.title = ggplot2::element_text(colour = "black"),
      axis.text = ggplot2::element_text(colour = "black"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", colour = "black"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = base - 1, colour = "black"),
      legend.text = ggplot2::element_text(size = base - 1, colour = "black"),
      legend.key = ggplot2::element_blank(),
      legend.key.size = ggplot2::unit(0.32, "cm"),
      legend.box.spacing = ggplot2::unit(1, "pt"),
      plot.title = ggplot2::element_text(face = "bold", colour = "black",
                                         hjust = 0.5, size = base + 1),
      plot.subtitle = ggplot2::element_text(colour = "grey35", hjust = 0.5,
                                            size = base - 1),
      panel.grid = ggplot2::element_blank()
    )
}

# 保存：PNG(300dpi) + 矢量 PDF + TIFF(LZW, 300dpi) + EPS(矢量)。文件前缀含期刊目录。
save_journal <- function(base, plot, width_mm, height_mm, dir = "figures/journal") {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  w <- width_mm / 25.4; h <- height_mm / 25.4
  fpng <- file.path(dir, paste0(base, ".png"))
  fpdf <- file.path(dir, paste0(base, ".pdf"))
  ftif <- file.path(dir, paste0(base, ".tiff"))
  feps <- file.path(dir, paste0(base, ".eps"))
  ggplot2::ggsave(fpng, plot, width = w, height = h, dpi = 300, bg = "white")
  ggplot2::ggsave(fpdf, plot, width = w, height = h, device = grDevices::cairo_pdf)
  ggplot2::ggsave(ftif, plot, width = w, height = h, dpi = 300, bg = "white",
                  device = ragg::agg_tiff, compression = "lzw")
  ggplot2::ggsave(feps, plot, width = w, height = h, device = grDevices::cairo_ps)
  cat("saved:", basename(fpng), basename(fpdf), basename(ftif), basename(feps), "\n")
  invisible(c(fpng, fpdf, ftif, feps))
}

# 色盲安全聚焦色板（与原图一致：TIER1 深红 #C8102E, TIER2 藏青 #005BBB, TIER3 灰）
tier_colours <- c(TIER1 = "#C8102E", TIER2 = "#005BBB", TIER3 = "#7F7F7F")
cohort_colours <- c(WIDE = "#1B9E77", INCLAVO = "#7570B3")
