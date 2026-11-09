# journal_fig5.R —— 论文主图 Fig5：PheWAS 多效性 + Steiger 方向性
# a) 5 候选 IV PheWAS（2,754 表型，GW 线，OLP 高亮，facet 5 基因）
# b) Steiger 方向检验散点（z_exp vs z_out，5 候选标记）
suppressMessages({library(ggplot2); library(data.table); library(ggrepel)})
source("analysis/scripts/journal_theme.R")

d <- fread("results/tables/phewas_olp_all.csv")
P_GW <- 5e-8
d <- d[!is.na(pval) & is.finite(-log10(pval))]   # 仅可用统计的表型（每 IV 308-415）
d[, mlogp := -log10(pval)]
d[, x := seq_len(.N), by = .(gene)]   # 文件顺序已按 category_index, phenocode 排

# plot_class：OLP / GW same-dir / GW opp-dir / not GW
# Direction is defined against THAT GENE'S OWN OLP effect, on the same allele basis
# (FinnGen ALT). An earlier version tested `beta < 0` for every gene, which assumes all
# five OLP effects are negative; they are not (GDF11 and SUOX are positive), so four of
# the five panels were classified the wrong way round.
olp_terms <- c("ORAL_LICHEN", "LICHEN_PLANUS")
d[, is_olp := grepl(paste(olp_terms, collapse = "|"), phenocode, ignore.case = TRUE)]
olp_ref <- d[is_olp == TRUE][order(pval)][, .SD[1], by = gene][, .(gene, beta_olp = beta)]
d <- merge(d, olp_ref, by = "gene", all.x = TRUE)
d[, dir_vs_olp := fifelse(is_olp, "olp",
                     fifelse(sign(beta) == sign(beta_olp), "same", "opposite"))]
d[is_olp == TRUE, plot_class := "OLP"]
d[is.na(plot_class) & mlogp >= -log10(P_GW) & dir_vs_olp == "opposite", plot_class := "GW opp-dir"]
d[is.na(plot_class) & mlogp >= -log10(P_GW) & dir_vs_olp == "same", plot_class := "GW same-dir"]
d[is.na(plot_class), plot_class := "not GW"]
d[, plot_class := factor(plot_class,
                         levels = c("OLP", "GW same-dir", "GW opp-dir", "not GW"))]

top3 <- d[plot_class %in% c("GW same-dir", "GW opp-dir")][order(pval)][, head(.SD, 3), by = .(gene)]
lab <- unique(rbind(d[is_olp == TRUE, .(gene, phenocode, pval, mlogp, plot_class, x)],
                    top3[, .(gene, phenocode, pval, mlogp, plot_class, x)]),
              by = c("gene", "phenocode"))
# 缩短长表型名以缓解拥挤分面内重叠
lab[, phenocode := sub("K11_ORAL_LICHEN_PLANUS_INCLAVO", "OLP_INCLAVO", phenocode)]
lab[, phenocode := sub("K11_LICHEN_PLANUS_WIDE", "LP_WIDE", phenocode)]
lab[, phenocode := sub("K11_ORAL_LICHEN_PLANUS", "OLP", phenocode)]

pal <- c(`OLP` = "#C8102E", `GW same-dir` = "#1B9E77", `GW opp-dir` = "#2B6CB0",
         `not GW` = "#B0B0B0")
shp <- c(`OLP` = 18, `GW same-dir` = 16, `GW opp-dir` = 17, `not GW` = 1)

pa <- ggplot(d, aes(x = x, y = mlogp)) +
  geom_hline(yintercept = -log10(P_GW), linetype = "dashed", colour = "grey45", linewidth = 0.3) +
  geom_point(aes(colour = plot_class, shape = plot_class), size = 0.55, alpha = 0.75) +
  geom_point(data = d[plot_class == "OLP"], aes(colour = plot_class, shape = plot_class), size = 1.1) +
  geom_text_repel(data = lab, aes(label = phenocode, colour = plot_class),
                  size = 1.7, max.overlaps = 12, seed = 1, min.segment.length = 0.2,
                  force = 0.8, box.padding = 0.25, max.time = 1.5) +
  facet_wrap(~ gene, ncol = 5, scales = "fixed") +
  scale_colour_manual(values = pal, name = NULL) +
  scale_shape_manual(values = shp, name = NULL) +
  scale_x_continuous(breaks = c(1, 200, 400), expand = expansion(mult = c(0.005, 0.01))) +
  scale_y_continuous(limits = c(0, 325), breaks = c(0, 100, 200, 300)) +
  labs(x = "FinnGen R13 phenotypes with available association (ordered by category)",
       y = expression(-log[10](italic(P))),
       title = "PheWAS: all 5 candidate IVs are pleiotropic across FinnGen R13 phenotypes",
       subtitle = paste0("dashed line = genome-wide 5e-8; red diamond = OLP endpoint; ",
                         "same/opposite are relative to that gene's own OLP effect; ",
                         "usable phenotypes per IV: 415/308/392/298/310 (facet order)")) +
  theme_journal(8) +
  theme(axis.text.x = element_text(size = 5.5),
        legend.text = element_text(size = 6.5),
        plot.margin = margin(6, 8, 6, 6))

# b) Steiger 方向性：z_exp vs z_out（eQTL vs OLP WIDE），5 候选标注
st <- fread("results/tables/steiger_direction_test.csv")
stwide <- st[grepl("LICHEN_PLANUS", pheno)]
stwide[, is_cand := fifelse(gene %in% c("GDF11","SUOX","EGFL8","TAPBP","HLA-DPA1"), gene, "")]
pb <- ggplot(stwide, aes(x = z_exp, y = z_out)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(aes(colour = is_cand == ""), size = 1.3, alpha = 0.8) +
  geom_point(data = stwide[is_cand != ""], aes(x = z_exp, y = z_out), shape = 18, size = 3, colour = "#C8102E") +
  geom_text_repel(data = stwide[is_cand != ""], aes(label = is_cand), size = 2.6, seed = 1, colour = "#C8102E") +
  scale_colour_manual(values = c(`TRUE` = "#9AA5B1", `FALSE` = "#C8102E"), guide = "none") +
  labs(x = "eQTLGen cis-eQTL Z (expression)", y = "FinnGen OLP WIDE MR Z (outcome)",
       title = "Steiger direction test: forward orientation for 4 of 5 candidates",
       subtitle = paste0("red = tiered candidate; EGFL8 forward but P = 0.082; ",
                         "GDF11 does not pass in INCLAVO (P = 0.72)")) +
  theme_journal(8)

suppressMessages(library(patchwork))
fig5 <- (pa + labs(tag = "A")) / (pb + labs(tag = "B")) &
  ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 11))
save_journal("fig5_phewas_steiger", fig5, 180, 170)
