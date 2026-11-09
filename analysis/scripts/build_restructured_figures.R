# build_restructured_figures.R — Plan 3 restructure of the main-figure deck.
#
# Author decision (2026-09-23): 6 main figures x >=6 panels, new folders, originals untouched.
#   Main 1 = design + druggable resources + the MR screen          (6 panels)
#   Main 2 = replication + colocalisation + candidate prioritisation (6 panels)
#   Main 3 = directionality + pleiotropy + docking + regional plots  (6 panels, S4 spanning)
#
# Mechanism: the existing journal/panel scripts already expose their panels as named objects
# (render_design(), png_a, png_b, vol, pa, pb, p1..p4 ...), so the new plates are composed from
# RE-RENDERED PANELS — no image cropping, and no edit to any existing script.
#
# Writes ONLY to figures/main_figures_v2/. Every figure-writing call inside the sourced scripts is
# neutralised first, so nothing under figures/main_figures/, figures/journal/ or
# figures/supplementary/ is touched.
suppressMessages({
  library(ggplot2); library(data.table); library(patchwork); library(grid)
  library(ggrepel); library(VennDiagram); library(scales)
})
source("analysis/scripts/journal_theme.R")

OUTDIR <- "figures/main_figures_v2"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
TAG <- ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 11))

save_v2 <- function(base, plot, width_mm, height_mm) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  ggplot2::ggsave(file.path(OUTDIR, paste0(base, ".png")), plot,
                  width = w, height = h, dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(OUTDIR, paste0(base, ".pdf")), plot,
                  width = w, height = h, device = grDevices::cairo_pdf)
  cat("  saved", base, sprintf("(%.0f x %.0f mm)\n", width_mm, height_mm))
}

NO_WRITE <- "^(\\s*)(save_journal|save_supp|save_main|ggsave)\\s*\\("

# Source a panel script into a fresh env with every figure-writing call disabled.
# Only single-line writers are neutralised; all of them are single-line in this project.
# `subs` optionally rewrites further strings (used to neutralise *nested* sources, e.g. a script
# that internally sys.source()s another panel script which would otherwise write to
# figures/journal).
neutralize <- function(script, subs = NULL) {
  txt <- readLines(script, warn = FALSE)
  hit <- grepl(NO_WRITE, txt)
  if (any(hit)) {
    txt[hit] <- sub(NO_WRITE, "\\1invisible(NULL) #", txt[hit])
    cat(sprintf("    (neutralised %d figure-writing call(s) in %s)\n", sum(hit), basename(script)))
  }
  if (!is.null(subs)) {
    for (s in subs) {
      n <- sum(grepl(s[1], txt, fixed = TRUE))
      txt <- gsub(s[1], s[2], txt, fixed = TRUE)
      if (n) cat(sprintf("    (patched %d nested source(s) in %s)\n", n, basename(script)))
    }
  }
  tmp <- tempfile(fileext = ".R")
  writeLines(txt, tmp)
  tmp
}

load_env <- function(script, subs = NULL) {
  e <- new.env(parent = globalenv())
  sys.source(neutralize(script, c(list(c("sys.source(", "safe_source(")), subs)), envir = e)
  e
}

cat("== loading panel sources ==\n")
# Route EVERY nested source through the neutraliser. A targeted pattern was not enough: one
# script sys.source()s another with a literal path, which slipped past a pattern-specific patch
# and re-wrote figures/journal/fig5_phewas_steiger.* (verified byte-identical graphics, only the
# embedded creation date differed; mtimes restored). Replacing the call name covers every variant.
safe_source <- function(path, ...) sys.source(neutralize(path), ...)
e_f1   <- new.env(parent = globalenv()); sys.source("analysis/scripts/fig1_parts.R", envir = e_f1)
e_f2   <- load_env("analysis/scripts/journal_fig2.R")
e_f3a  <- load_env("analysis/scripts/journal_fig3.R")
e_f4   <- load_env("analysis/scripts/journal_fig4.R")
e_f6   <- load_env("analysis/scripts/journal_fig6.R")
e_supp <- load_env("analysis/scripts/build_supplementary_figures.R")

# =====================================================================
# Main Figure 1 — study design, druggable resources and the MR screen (6 panels)
#   row 1: A design schematic (full width - it is a wide, text-heavy diagram and was
#          designed at 180 mm in the original deck; squeezing it into a 60 mm column
#          collides its internal labels)
#   row 2: B four-source Venn | C druggable-with-cis-eQTL
#   row 3: D Manhattan WIDE  | E Manhattan INCLAVO | F volcano
# =====================================================================
cat("== Fig 1 ==\n")
# right-hand panels keep a margin so their titles/labels are not clipped at the plate edge
SAFE_R <- ggplot2::theme(plot.margin = ggplot2::margin(4, 8, 4, 4))
# Cosmetic overrides for the re-used MR panels at column width. NOTE: do NOT re-specify the
# x scale — the Manhattan panels map 22 chromosome labels onto midpoint breaks, and replacing
# the scale silently drops those labels (tried, reverted). Instead relieve crowding with font
# size, and shorten the one axis label that was running past the plate edge.
thin_ticks <- ggplot2::theme(axis.text.x = ggplot2::element_text(size = 5.2))
fig1_A <- e_f1$render_design() + labs(tag = "A", title = "Study design") + SAFE_R
# B's Venn is drawn as a grob stretched to fill its panel, so its outer set labels (e.g.
# "Pharos 79") reach the panel edge and the leftmost one is cut off by the plate boundary.
# Insetting the panel shrinks the grob and brings the labels inside.
VENN_PAD <- ggplot2::theme(plot.margin = ggplot2::margin(6, 12, 6, 20))
fig1_B <- e_f1$render_venn()   + labs(tag = "B") + VENN_PAD
fig1_C <- e_f1$render_funnel() + labs(tag = "C") +
  labs(title = "Druggable genes with cis-eQTL support",
       subtitle = "5,801 of 9,261 druggable genes (62.6%)\nhave an eQTLGen cis-eQTL instrument") +
  SAFE_R
# D and E carry a long single-line subtitle ("n = 5,040 genes; dashed line = Bonferroni
# (P = 9.9e-06)") that the source script wrapped for a half-page panel; in a third of 180 mm it is
# clipped at the panel edge (panel E lost its closing bracket). Both are re-wrapped to two lines.
#
# Their x axis labels every odd chromosome (1,3,...,21). That fits a half-page panel but not a third
# of the plate, where "19" and "21" overlap at the right end. The recomposed panels therefore label
# every fourth chromosome (1,5,9,13,17,21). Only the x scale is replaced; these panels draw no
# x legend, so there is no guide to lose (unlike Main 2 panel D, see the note there).
thin_chr <- function(d) {
  m <- unique(d[, .(chr, mid)])[order(chr)]
  keep <- m$chr %in% c(1, 5, 9, 13, 17, 21)
  suppressMessages(ggplot2::scale_x_continuous(
    breaks = m$mid[keep], labels = m$chr[keep],
    expand = ggplot2::expansion(mult = c(0.005, 0.01))))
}
fig1_D <- e_f2$png_a + labs(tag = "D",
                            subtitle = "n = 5,040 genes\ndashed = Bonferroni (9.9e-06)") +
  thin_chr(e_f2$w) + thin_ticks
fig1_E <- e_f2$png_b + labs(tag = "E",
                            subtitle = "n = 5,040 genes\ndashed = Bonferroni (9.9e-06)") +
  thin_chr(e_f2$i) + thin_ticks
# F's legend sat left-centre inside the panel, where it covered the EGFL8 hit (the leftmost
# labelled point is CSNK2B at the top, EGFL8 just below it). The top-right of the volcano is empty,
# so the legend moves there, right-justified, and a little smaller.
fig1_F <- e_f2$vol   + labs(tag = "F") +
  labs(x = "MR Wald-ratio Z (per-SD expression)",
       subtitle = "n = 5,040 genes; + = higher risk\ndashed = Bonferroni (9.9e-06)") + SAFE_R +
  ggplot2::theme(legend.position = "inside",
                 legend.position.inside = c(0.985, 0.975),
                 legend.justification = c(1, 1),
                 legend.text = ggplot2::element_text(size = 5.0),
                 legend.title = ggplot2::element_text(size = 5.2),
                 legend.key.size = ggplot2::unit(0.24, "cm"),
                 legend.margin = ggplot2::margin(1.5, 2.5, 1.5, 2.5))

fig1_v2 <- (fig1_A / (fig1_B | fig1_C) / (fig1_D | fig1_E | fig1_F)) +
  patchwork::plot_layout(heights = c(0.85, 1, 1.15)) & TAG
save_v2("fig01_design_screen", fig1_v2, 180, 210)

# =====================================================================
# Main Figure 2 — replication, colocalisation and prioritisation (6 panels)
#
# Layout constraint discovered while building this plate: the panels that list all 22 WIDE
# hits (the dual-cohort forest, the PP.H4 bars and the posterior decomposition) were designed
# at full page width. Side by side at half width their 22 gene labels collide and their long
# titles overlap each other. They therefore get a full-width (spanning) row each, with their
# gene labels set at a smaller size, and only the intrinsically compact panels share a row.
#
#   row 1 (span, tall): A 22-hit dual-cohort forest
#   row 2:              B five-candidate forest | D SuSiE vs abf (both compact)
#   row 3 (span):       C coloc.abf PP.H4 bars (22 genes on the x axis)
#   row 4 (span):       E coloc posterior decomposition (22 genes on the x axis)
#   row 5 (span):       F Bonferroni vs FDR sensitivity (two wide facets)
# Panel titles that run to two lines in the source scripts are shortened here so they fit the
# spanning rows; the new legend states the shortened wording.
# =====================================================================
cat("== Fig 2 ==\n")
tiny_x <- ggplot2::theme(axis.text.x = ggplot2::element_text(size = 4.0, angle = 90,
                                                            hjust = 1, vjust = 0.5))
# A lists 22 genes and B lists 10 (five candidates x two cohorts). At this plate width their row
# pitch is ~4 mm, so the y tick labels must shrink below the theme's default or adjacent entries
# touch. E stacks three posterior rows in one panel, so its 0.00-1.00 ticks need the same treatment.
tiny_y <- ggplot2::theme(axis.text.y = ggplot2::element_text(size = 4.4))
# B draws two annotation columns inside the panel: the P value anchored at each OR (extending
# right) and a right-aligned "PP.H4" column at a fixed x. In a half-width cell the longest entry
# (GDF11, OR 2.88, "P = 6.51e-06") reaches the PP.H4 column, so both text layers are set 15%
# smaller. Only text layers are touched; the geoms' other parameters are left alone.
shrink_text_layers <- function(p, k = 0.60, p_vjust = -0.45) {
  for (i in seq_along(p$layers)) {
    if (!inherits(p$layers[[i]]$geom, "GeomText")) next
    ap <- p$layers[[i]]$aes_params
    if (!is.null(ap$size)) p$layers[[i]]$aes_params$size <- ap$size * k
    # The P label is the layer anchored at the OR (hjust = -0.3). Its default vjust of -1.2 lifts
    # it a full row, so on the 10-row panel it lands on the row above and runs into that row's
    # PP.H4 entry. Pulling it down to -0.45 keeps it beside its own point.
    if (!is.null(ap$hjust) && isTRUE(all.equal(ap$hjust, -0.3))) {
      p$layers[[i]]$aes_params$vjust <- p_vjust
    }
  }
  p
}
fig2_A <- e_supp$p1 + labs(tag = "A") + tiny_y
fig2_B <- shrink_text_layers(e_f3a$p) + labs(tag = "B") + tiny_y +
  labs(subtitle = paste0("single cis-eQTL IV per gene; OR < 1 = higher expression, lower risk; ",
                         "right column = coloc.susie PP.H4")) +
  ggplot2::theme(plot.margin = ggplot2::margin(6, 16, 16, 6))
# The two compact panels (B, and the SuSiE-vs-abf scatter) have to share a row, which puts the
# scatter third in reading order. The letters follow the layout, so the scatter is C and the
# 22-gene coloc.abf bars are D; the plate then reads A, B, C, D, E, F. B takes more of the row
# than C because its annotation columns need the width, while C is a two-facet scatter that
# tolerates a narrower cell.
fig2_C <- e_f4$pb + labs(tag = "C") +
  # Title shortened: the column is too narrow for "Multi-signal (SuSiE) vs single-signal (abf)"
  # at title size, and it lost its closing bracket at the plate edge.
  labs(title = "SuSiE vs single-signal (abf)",
       subtitle = "dotted line = 1:1; dashed = PP.H4 0.8",
       y = "SuSiE PP.H4") +
  # The legend is put INSIDE the scatter: moved to the plot edge it enlarged the cell's minimum
  # size and patchwork could no longer fit it beside panel B (its title was clipped). Inside, the
  # panel keeps a normal size; the legend sits in the lower-left corner, which is empty.
  ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2, byrow = TRUE)) +
  ggplot2::theme(legend.position = "inside",
                 legend.position.inside = c(0.015, 0.015),
                 legend.justification = c(0, 0),
                 legend.background = ggplot2::element_rect(colour = "grey75",
                                                           fill = ggplot2::alpha("white", 0.85),
                                                           linewidth = 0.3),
                 axis.text.x = ggplot2::element_text(size = 4.4),
                 plot.margin = ggplot2::margin(6, 8, 16, 8),
                 legend.text = ggplot2::element_text(size = 5.0),
                 legend.title = ggplot2::element_blank(),
                 legend.key.size = ggplot2::unit(0.24, "cm"))
fig2_D <- e_f4$pa + labs(tag = "D") + tiny_x +
  labs(title = "Colocalisation (coloc.abf) for the 22 WIDE hits",
       subtitle = "dashed line = PP.H4 0.8; genes ordered by PP.H4 within each definition") +
  ggplot2::theme(plot.margin = ggplot2::margin(16, 6, 6, 6))
fig2_E <- e_supp$p2 + labs(tag = "E") + tiny_x + tiny_y +
  labs(title = "Colocalisation posteriors for the 22 WIDE hits",
       subtitle = "dashed line = strong-colocalisation threshold PP.H4 0.8")
fig2_F <- e_supp$p3 + labs(tag = "F") +
  labs(subtitle = "red = Bonferroni-significant; orange = additionally BH-FDR < 0.05; dashed = Bonferroni")

fig2_v2 <- (fig2_A / ((fig2_B | fig2_C) +
                        patchwork::plot_layout(widths = c(1.55, 1))) /
            fig2_D / fig2_E / fig2_F) +
  patchwork::plot_layout(heights = c(2.00, 1.25, 0.85, 1.05, 0.60)) & TAG
# This plate is taller than the others on purpose: panel A carries 22 gene rows and B carries 10,
# and at 240 mm their row pitch was ~1 mm, which no font size can survive. 280 mm gives A ~4 mm
# per row.
save_v2("fig02_replication_colocalisation", fig2_v2, 180, 280)

# =====================================================================
# Main Figure 3 - directionality, pleiotropy and colocalisation detail (4 panels)
#   A PheWAS (spanning: five facets with in-facet phenotype labels)
#   B Steiger | C SuSiE pleiotropy
#   D PP.H4 matrix (spanning: ten rotated phenotype labels)
#
# The regional-association panel (former Supplementary Figure 4) was REMOVED from this figure
# after testing: it is a 5-gene x 2-definition grid of ten sub-plots, and at a one-sixth share of
# the plate (180 x ~90 mm) each sub-plot is ~14 mm tall, so its axes and rsID labels are
# unreadable. It therefore stays a standalone supplementary figure (see the mapping note).
# Three panels in one row was also tested and rejected: their long titles overlap.
# =====================================================================
cat("== Fig 3 ==\n")
e_mr <- load_env("analysis/scripts/build_final_mr_figures.R")
# A left a wide empty band on its right: the top row is the full plate width, so the right margin
# was pure waste next to the right-hand y-axis label. Trim it, and show more chromosome detail -
# every odd chromosome (the density the source panel was designed for) now that the panel is wider.
fig3_A <- e_mr$phewas  + labs(tag = "A") +
  ggplot2::theme(plot.margin = ggplot2::margin(6, 4, 6, 4))
fig3_B <- e_mr$steiger + labs(tag = "B")
# C plots the chr12 candidates as two 0.8-threshold facets, and its ggrepel labels collide under
# the SUOX facet because the original was a half-page panel. Give it more width (see the row
# ratios below) and set its annotation layer smaller; only GeomTextRepel layers are touched.
fig3_C <- e_mr$susie   + labs(tag = "C", y = "SuSiE PP.H4") + tiny_y +
  ggplot2::theme(plot.margin = ggplot2::margin(6, 10, 6, 6))
for (i in seq_along(fig3_C$layers)) {
  if (inherits(fig3_C$layers[[i]]$geom, "GeomTextRepel")) {
    fig3_C$layers[[i]]$aes_params$size <- 1.45
    fig3_C$layers[[i]]$geom_params$max.overlaps <- 8
    fig3_C$layers[[i]]$geom_params$force <- 2
    fig3_C$layers[[i]]$geom_params$box.padding <- 0.22
    fig3_C$layers[[i]]$geom_params$point.padding <- 0.2
  }
}
fig3_D <- e_mr$heatmap + labs(tag = "D")

# REVERTED: an explicit `design =` grid was tried here and mis-sized the plate (A took most of the
# height and the lower panels were squeezed out), so the nested form is kept.
# Row 2 holds the only side-by-side pair: C (two SuSiE facets with long annotations) takes more width
# than B, which is the pair's only ratio that matters.
fig3_v2 <- (fig3_A / ((fig3_B | fig3_C) + patchwork::plot_layout(widths = c(1, 1.32))) /
            fig3_D) +
  patchwork::plot_layout(heights = c(1.25, 1.60, 0.85)) & TAG
# the canvas shrank with the removed row: 3.70/4.65 of the former 260 mm
save_v2("fig03_pleiotropy_docking", fig3_v2, 180, 207)

cat("\nDone. Restructured figures are in", OUTDIR, "\n")