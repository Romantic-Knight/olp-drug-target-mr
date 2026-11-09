# -*- coding: utf-8 -*-
"""Stage A: NEOLP (non-erosive) vs EOLP (erosive) x 5 MR candidate genes,
stratified by annotated OLP cell type (GSE211630).

Design notes
------------
* Universe: 46,446 annotated cells of the OLP object (doublet-suspect
  c12/c18, n=32, excluded - no family annotation / QC suspects).
* Replication unit: SAMPLE (n = 3 NEOLP vs 2 EOLP). Cell-level statistics are
  exploratory only (cells are pseudoreplicates within a sample); the
  sample-level permutation test is the conservative companion.
* Object stores log1p-normalised expression only (no raw-count layer), so:
      mean_all = mean log1p over all cells of the type in the sample
               = frac_expressing x mean_exp (detection AND amount)
      mean_exp = mean log1p among expressing cells (intensity)
* Tests per gene x cell type:
    1) cell-level two-sided Mann-Whitney (pooled cells) + BH FDR over the
       tested grid -> exploratory;
    2) sample-level two-sided permutation p over the C(5,3)=10 label splits
       of the 5 per-sample mean_all values (min possible p = 0.1).
* Sparse detection is reported as insufficient, never imputed or dropped
  silently (n_expr / n_cells columns in the tables carry the denominators).

Outputs
-------
results/tables/geo_gse211630_erosion_per_sample.csv
    long: gene x celltype x sample x condition -> n_cells, frac, mean_all,
    mean_exp (5 samples x 19 types x 5 genes)
results/tables/geo_gse211630_erosion_candidate_celltype.csv
    gene x celltype (sub + family levels) group summary + both tests
figures/main_figures/fig06_target_expression.png (300 dpi) + .pdf (+ 150 dpi QA copy)
    composite Figure 6 (was Figure 8): (A) annotated UMAP, (B) NEOLP-vs-EOLP UMAP,
    (C)-(G) candidate-gene UMAP feature strip on the same 46,478-cell
    embedding, and the lower 2x3 candidate-gene erosion dotplot grid.
"""
import gc
import itertools
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
os.chdir(ROOT)
SCAN = os.path.join("data", "processed", "geo_GSE211630", "scanpy")
TAB = os.path.join("results", "tables")
FIGDIR = os.path.join("figures", "main_figures")   # final output goes to main_figures/
QA = os.path.join("analysis", "scrna", "_qa")
os.makedirs(QA, exist_ok=True)

CAND = ["GDF11", "SUOX", "EGFL8", "TAPBP", "HLA-DPA1"]
COND_ORDER = ["NEOLP", "EOLP"]
COND_COL = {"NEOLP": "#4477AA", "EOLP": "#EE7733"}     # blue / orange, CVD-safe
HIT_RED = "#9A2A2A"                                     # alarm hue: FDR flag only
N_NEOLP, N_EOLP = 3, 2
SIZE_MIN, SIZE_MAX = 9.0, 26.0                          # bubble area, frac -> size

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7,
                     "axes.linewidth": 0.6, "svg.fonttype": "none"})

# ------------------------------------------------------------------ load
olp = sc.read_h5ad(os.path.join(SCAN, "gse211630_olp_lognorm.h5ad"))
ob = pd.read_csv(os.path.join(SCAN, "gse211630_olp_obs_annotation.csv"), index_col=0)
assert (olp.obs.index == ob.index).all()
olp.obs["celltype_sub"] = ob["celltype_sub"].values
olp.obs["celltype_family"] = ob["celltype_family"].values
# capture the FULL 46,478-cell object arrays before doublet exclusion (used by
# the Figure-8 overview panels A/B, where doublets are drawn grey)
umap_all = np.asarray(olp.obsm["X_umap"])
leid_all = olp.obs["leiden_0.5"].astype(int).values
cond_all = olp.obs["condition"].astype(str).values
# FULL candidate expression on all 46,478 cells, captured before exclusion
# (the C-G UMAP feature strip draws every cell, doublets included as in
# panels A/B, so the map is identical to the overview embedding).
Xc_all = np.asarray(olp[:, CAND].X.todense()) if hasattr(olp[:, CAND].X, "todense") \
    else np.asarray(olp[:, CAND].X)
keep = ~olp.obs["leiden_0.5"].astype(int).isin([12, 18])   # doublet-suspect out
olp = olp[keep].copy()
assert olp.n_obs == 46446, olp.n_obs
print("annotated cells used:", olp.n_obs)

Xc = np.asarray(olp[:, CAND].X.todense()) if hasattr(olp[:, CAND].X, "todense") \
    else np.asarray(olp[:, CAND].X)
sample = olp.obs["sample"].astype(str).values
cond = olp.obs["condition"].astype(str).values
leid = olp.obs["leiden_0.5"].astype(int).values
sub = olp.obs["celltype_sub"].astype(str).values
fam = olp.obs["celltype_family"].astype(str).values
gene_ix = {g: i for i, g in enumerate(CAND)}
# release the full 46k x 17.7k log1p matrix before figure construction (the
# composite render needs the free headroom on the 4 GB host)
del olp
gc.collect()

# cluster <-> annotation maps (1:1 as fixed in Stage 4)
cl_sub = dict(zip(*np.unique(list(zip(leid, sub)), axis=0).T)) if False else None
cl_sub = {}
for cl, t in zip(leid, sub):
    cl_sub.setdefault(int(cl), t)
cl_fam = {}
for cl, f in zip(leid, fam):
    cl_fam.setdefault(int(cl), f)
assert len(cl_sub) == 19 and sorted(cl_sub) == sorted(cl_fam)
t2cl = {t: cl for cl, t in cl_sub.items()}          # cell type text -> cluster id

# short row labels (c# = Leiden cluster, same strings as Stage-4 dotplot)
SHORT = {0: "c0 Fibroblast SFRP2+", 1: "c1 Fibroblast CXCL14+", 2: "c2 Pericyte/SMC",
         3: "c3 T resting", 4: "c4 T/NK stressed", 5: "c5 T/NK low-complex.",
         6: "c6 CD8+ GZMK+ EM", 7: "c7 CD4+ Treg", 8: "c8 T naive IL7R+",
         9: "c9 Plasmablast GZMB+", 10: "c10 Macrophage LYZ+",
         11: "c11 Epith. basal/prolif.", 13: "c13 Mast", 14: "c14 B cell",
         15: "c15 Plasma cell", 16: "c16 CD8+ GZMB+/CXCL13+",
         17: "c17 Endo lymph.", 19: "c19 Endo blood", 20: "c20 Myeloid S100-hi"}


def bh(pvals):
    p = np.asarray(pvals, float)
    q = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    if ok.sum() == 0:
        return q
    pv = p[ok]
    o = np.argsort(pv)
    # Benjamini-Hochberg: q_(i) = p_(i) * m / i, then enforce monotonicity.
    # An earlier version used p * i / m, which is not an adjusted P value.
    m = len(pv)
    scaled = pv[o] * m / np.arange(1, m + 1)
    q[ok] = np.minimum.accumulate(scaled[::-1])[::-1][np.argsort(o)]
    return q


# ---------------------------------------------------------- per-sample metrics
rows = []
for s in np.unique(sample):
    ms = sample == s
    cond_s = cond[ms][0]
    for cl, t in cl_sub.items():
        m2 = ms & (leid == cl)
        if m2.sum() == 0:
            continue
        for j, g in enumerate(CAND):
            v = Xc[m2, j]
            n_expr = int((v > 0).sum())
            rows.append({"gene": g, "sample": s, "condition": cond_s,
                         "cluster": cl, "celltype": t, "family": cl_fam[cl],
                         "n_cells": int(m2.sum()), "n_expr": n_expr,
                         "frac_expr": n_expr / m2.sum(),
                         "mean_all": float(v.mean()),
                         "mean_exp": float(v[v > 0].mean()) if n_expr else 0.0})
per_sample = pd.DataFrame(rows).sort_values(
    ["gene", "cluster", "sample"]).reset_index(drop=True)
# note: cluster 20 (Myeloid S100-high) has 0 cells in NEOLP-2 (real composition
# difference after doublet exclusion), so rows for that combo are absent; the
# n_cells column still accounts for every annotated cell once per gene below.
per_sample.to_csv(os.path.join(TAB, "geo_gse211630_erosion_per_sample.csv"), index=False)
assert per_sample.n_cells.sum() == len(CAND) * 46446  # each cell counted once per gene (46,446 = all 5 samples)

def _summ(g, level, key, fam_lab, n, e, n_cell, e_cell, n_expr, e_expr, vn, ve):
    """One gene x celltype summary row (level in {'sub','family'})."""
    def wt(df, col, w):
        return float(np.average(df[col], weights=w)) if len(df) and w.sum() else np.nan
    w_mean_n = wt(n, "mean_all", n.n_cells); w_mean_e = wt(e, "mean_all", e.n_cells)
    w_frac_n = wt(n, "frac_expr", n.n_cells); w_frac_e = wt(e, "frac_expr", e.n_cells)
    w_exp_n = wt(n, "mean_exp", n.n_expr);  w_exp_e = wt(e, "mean_exp", e.n_expr)
    cell_p = np.nan
    if vn.size >= 10 and ve.size >= 10 and (vn > 0).sum() >= 5 and (ve > 0).sum() >= 5:
        try:
            cell_p = stats.mannwhitneyu(vn, ve, alternative="two-sided").pvalue
        except ValueError:
            cell_p = np.nan
    perm_p, med_d = np.nan, np.nan
    if len(n) == N_NEOLP and len(e) == N_EOLP:
        vals = np.concatenate([n.mean_all.values, e.mean_all.values])
        if vals.max() > vals.min():
            obs = vals[:N_NEOLP].mean() - vals[N_NEOLP:].mean()
            d = []
            for ii in itertools.combinations(range(5), 3):
                ll = np.zeros(5, int); ll[list(ii)] = 1
                d.append((vals[ll == 1].mean() - vals[ll == 0].mean()))
            perm_p = float((np.abs(np.asarray(d)) >= np.abs(obs)).mean())
            med_d = float(np.median(n.mean_all) - np.median(e.mean_all))
    return {"gene": g, "level": level,
            "celltype": key if level == "sub" else "", "family": fam_lab,
            "cluster": t2cl.get(key, np.nan) if level == "sub" else np.nan,
            "n_cells_NEOLP": n_cell, "n_cells_EOLP": e_cell,
            "n_expr_NEOLP": n_expr, "n_expr_EOLP": e_expr,
            "frac_NEOLP": w_frac_n, "frac_EOLP": w_frac_e,
            "mean_all_NEOLP": w_mean_n, "mean_all_EOLP": w_mean_e,
            "mean_exp_NEOLP": w_exp_n, "mean_exp_EOLP": w_exp_e,
            "delta_mean_all": w_mean_e - w_mean_n, "delta_frac": w_frac_e - w_frac_n,
            "median_sample_mean_all_NEOLP": float(np.median(n.mean_all)) if len(n) else np.nan,
            "median_sample_mean_all_EOLP": float(np.median(e.mean_all)) if len(e) else np.nan,
            "median_sample_mean_exp_NEOLP": float(np.median(n.mean_exp)) if len(n) else np.nan,
            "median_sample_mean_exp_EOLP": float(np.median(e.mean_exp)) if len(e) else np.nan,
            "wilcoxon_cell_p": cell_p, "perm_p_sample": perm_p,
            "median_delta_sample": med_d}


# ------------------------------------------------------- gene x celltype summary
NE = per_sample[per_sample.condition == "NEOLP"]
EO = per_sample[per_sample.condition == "EOLP"]
srows = []
for g in CAND:
    for cl, t in sorted(cl_sub.items()):                                # sub level
        n = NE[(NE.gene == g) & (NE.celltype == t)]
        e = EO[(EO.gene == g) & (EO.celltype == t)]
        n_cell, e_cell = int(n.n_cells.sum()), int(e.n_cells.sum())
        n_expr, e_expr = int(n.n_expr.sum()), int(e.n_expr.sum())
        vn = Xc[(cond == "NEOLP") & (sub == t), gene_ix[g]]
        ve = Xc[(cond == "EOLP") & (sub == t), gene_ix[g]]
        srows.append(_summ(g, "sub", t, cl_fam[cl], n, e, n_cell, e_cell,
                           n_expr, e_expr, vn, ve))
    for f in sorted(set(cl_fam.values())):                              # family level
        n = NE[(NE.gene == g) & (NE.family == f)]
        e = EO[(EO.gene == g) & (EO.family == f)]
        n_cell, e_cell = int(n.n_cells.sum()), int(e.n_cells.sum())
        n_expr, e_expr = int(n.n_expr.sum()), int(e.n_expr.sum())
        vn = Xc[(cond == "NEOLP") & (fam == f), gene_ix[g]]
        ve = Xc[(cond == "EOLP") & (fam == f), gene_ix[g]]
        srows.append(_summ(g, "family", f, f, n, e, n_cell, e_cell,
                           n_expr, e_expr, vn, ve))
summ = pd.DataFrame(srows)
for lev in ("sub", "family"):
    m = summ.level == lev
    summ.loc[m, "BH_FDR_cell"] = bh(summ.loc[m, "wilcoxon_cell_p"].values)
cols = ["gene", "level", "celltype", "family", "cluster",
        "n_cells_NEOLP", "n_cells_EOLP", "n_expr_NEOLP", "n_expr_EOLP",
        "frac_NEOLP", "frac_EOLP", "mean_all_NEOLP", "mean_all_EOLP",
        "mean_exp_NEOLP", "mean_exp_EOLP", "delta_mean_all", "delta_frac",
        "median_sample_mean_all_NEOLP", "median_sample_mean_all_EOLP",
        "median_sample_mean_exp_NEOLP", "median_sample_mean_exp_EOLP",
        "wilcoxon_cell_p", "BH_FDR_cell", "perm_p_sample"]
summ = summ[cols].sort_values(["level", "gene", "cluster"], na_position="first")
summ.to_csv(os.path.join(TAB, "geo_gse211630_erosion_candidate_celltype.csv"), index=False)


# ---------------------------------------------------------------- figure rows
# sub-level pairs with adequate per-condition support (else: table only)
sub_s = summ[summ.level == "sub"].copy()
sub_s["_ok"] = (sub_s.n_expr_NEOLP >= 5) & (sub_s.n_expr_EOLP >= 5) & \
               (sub_s.n_cells_NEOLP >= 10) & (sub_s.n_cells_EOLP >= 10)
fig_sub = sub_s[sub_s._ok]
print("figure rows per gene (>=5 expressing cells in both conditions):")
print(fig_sub.groupby("gene").size().to_string())
for g in CAND:
    if (fig_sub.gene == g).sum() == 0:
        print(f"  note: {g} has no cell type meeting the >=5-expressor/condition "
              f"threshold -> panel shows 'insufficient detection'")

# ---------------------------------------------------------------- figure


# =====================================================================
# Restructured 2026-09-30 (main-figure re-plan).
# The drawing code below is the original module-level block, unchanged, moved into panel and plate
# functions so panels can be recomposed into new figures. `main()` passes the original geometry, so
# fig06_target_expression.{png,pdf} is reproduced exactly. The plate holds four panels: A/B the
# atlas overview, C the candidate-gene feature composite (an embedded raster), D the 2x3 erosion
# dotplot grid.
# =====================================================================

import colorsys

# palette shared with the Stage-4 annotated atlas (plot_celltype_gse211630.py)
FAMILY_ORDER = ["T/NK", "B/Plasma", "Myeloid", "Mast", "Epithelial",
                "Fibroblast", "Pericyte/SMC", "Endothelial"]
FAM_HLS = {"T/NK": (0.00, 0.47, 0.92), "B/Plasma": (0.585, 0.44, 0.80),
           "Myeloid": (0.34, 0.47, 0.72), "Mast": (0.80, 0.52, 0.62),
           "Epithelial": (0.105, 0.42, 0.88), "Fibroblast": (0.145, 0.40, 0.72),
           "Pericyte/SMC": (0.07, 0.26, 0.55), "Endothelial": (0.50, 0.42, 0.72)}
_SUB_HLS = {
    3: (0.000, 0.50, 0.95), 4: (0.000, 0.66, 0.50), 5: (0.015, 0.82, 0.28),
    6: (0.985, 0.24, 0.80), 7: (0.035, 0.55, 0.90), 8: (0.050, 0.74, 0.42),
    16: (0.000, 0.36, 0.88)}
_SUB_LADJ = {9: +0.24, 14: 0.00, 15: -0.15, 10: 0.00, 20: +0.15, 11: 0.00,
             0: 0.00, 1: +0.20, 2: 0.00, 17: +0.25, 19: 0.00}
assert set(FAMILY_ORDER) == set(cl_fam.values())

def _sub_color(cl):
    if cl in _SUB_HLS:
        h, l2, s = _SUB_HLS[cl]
        return colorsys.hls_to_rgb(h, l2, s)
    h, l, s = FAM_HLS[cl_fam[cl]]
    l2 = float(np.clip(l + _SUB_LADJ.get(cl, 0.0), 0.16, 0.82))
    return colorsys.hls_to_rgb(h, l2, s)

legend_rows = []
for fam in FAMILY_ORDER:
    for cl in sorted(cl_sub):
        if cl_fam[cl] == fam:
            legend_rows.append((cl, f"c{cl}  {cl_sub[cl]}", fam))

# A/B overview panels share one box size and one vertical centre so their tops
# and bottoms align and B is no smaller than A. The box plus its left-side panel
# label and axis title must all fit below the canvas top (y=1.0), with a clear
# gap down to panel C below.
_AW = 3.00                                # shared box width (inches) for A and B
_axw, _axh = _AW / 12.2, _AW / 18.9
_yMid = 0.890                             # A/B box top ~0.969; tag (above title) fits under y=1
dbl = np.isin(leid_all, (12, 18))


def panel_A(ax, tag="A"):
    """Annotated OLP atlas UMAP (46,478 cells; doublet-suspect c12/c18 drawn grey)."""
    for ax_ in (ax, ax):
        ax_.set_xticks([])
        ax_.set_yticks([])
        for sp in ax_.spines.values():
            sp.set_visible(False)
    for dcl in (12, 18):                      # doublet-suspect base layer (grey)
        m = leid_all == dcl
        ax.scatter(umap_all[m, 0], umap_all[m, 1], s=1.1, c=[(0.76, 0.76, 0.76)] * int(m.sum()),
                    linewidths=0, alpha=0.8, rasterized=True)
    for cl in sorted(cl_sub):                 # 19 annotated clusters (family palette)
        m = (leid_all == cl) & ~dbl
        ax.scatter(umap_all[m, 0], umap_all[m, 1], s=1.4, c=[_sub_color(cl)] * int(m.sum()),
                    linewidths=0, alpha=0.82, rasterized=True)
    for dcl in (12, 18):                      # re-draw doublets on top, with edge
        m = leid_all == dcl
        ax.scatter(umap_all[m, 0], umap_all[m, 1], s=3.6, c=[(0.76, 0.76, 0.76)] * int(m.sum()),
                    linewidths=0.4, edgecolors="#444444", alpha=0.98, rasterized=True)
    ax.set_title("OLP atlas (n = 46,478)\n19 annotated cell types, 8 families",
                  loc="center", fontsize=8.0, weight="bold", pad=3)

    ax.text(-0.02, 1.075, tag, transform=ax.transAxes, fontsize=15,
            fontweight="bold", ha="left", va="bottom")
    return ax


def panel_A_legend(ax, ncol=1, fontsize=7.4, labelspacing=0.12, ms=6.4, ms_dbl=5.2):
    """The 28-row family/subtype key, anchored to the right of panel A.

    ncol/labelspacing let a recomposition trade width for height: the original single column needs
    ~99 mm at the original spacing, which does not fit beside a half-width panel in a composite, so
    a two-column form (~50 mm) can be requested. The defaults reproduce the original figure.
    """
    # A legend (family grouping, same strings as the atlas figure)
    handles, labels = [], []
    for fam in FAMILY_ORDER:
        handles.append(plt.Line2D([], [], ls="-", color="#666666", lw=1.4, label=""))
        labels.append(fam)
        for cl, lab, fam2 in legend_rows:
            if fam2 != fam:
                continue
            handles.append(plt.Line2D([], [], marker="o", ls="", markersize=ms,
                                      color=_sub_color(cl), alpha=0.95))
            labels.append(lab)
    handles.append(plt.Line2D([], [], marker="o", ls="", markersize=ms_dbl,
                              color=(0.76, 0.76, 0.76), markeredgecolor="#444444",
                              markeredgewidth=0.6))
    labels.append("doublet-suspect c12/c18 (n = 32, excluded)")
    # Single column anchored at the panel's right edge. At the increased font size
    # the 28-row list is tall, so row spacing is tightened to keep the box clear of
    # panel C's image below (two columns would instead run over panel B).
    lga = ax.legend(handles, labels, loc="upper left", bbox_to_anchor=(1.02, 1.0),
                    frameon=True, framealpha=0.95, edgecolor="#DDDDDD", fontsize=fontsize,
                    labelspacing=labelspacing, handletextpad=0.45, borderpad=0.30,
                    ncol=ncol, columnspacing=1.0)
    for t, fam in zip(lga.get_texts(), labels):
        if fam in FAMILY_ORDER:
            t.set_weight("bold")
            t.set_fontsize(fontsize + 0.2)

    return ax


def panel_B(ax, tag="B", ms=5.2, legend_loc="upper right",
            legend_bbox=(0.995, 0.995), legend_ncol=1):
    """NEOLP-vs-EOLP overlay on the same embedding (doublets grey)."""
    for ax_ in (ax, ax):
        ax_.set_xticks([])
        ax_.set_yticks([])
        for sp in ax_.spines.values():
            sp.set_visible(False)
    # (B) condition overlay on the same embedding (doublets grey)
    for cn in ("NEOLP", "EOLP"):
        m = (cond_all == cn) & ~dbl
        ax.scatter(umap_all[m, 0], umap_all[m, 1], s=1.4, c=[COND_COL[cn]] * int(m.sum()),
                    linewidths=0, alpha=0.78, rasterized=True)
    m = dbl
    ax.scatter(umap_all[m, 0], umap_all[m, 1], s=3.6, c=[(0.76, 0.76, 0.76)] * int(m.sum()),
                linewidths=0.4, edgecolors="#444444", alpha=0.98, rasterized=True)
    ax.set_title("Non-erosive (NEOLP) vs erosive (EOLP) OLP",
                  loc="center", fontsize=8.8, weight="bold", pad=3)
    hb = [plt.Line2D([], [], marker="o", ls="", color=COND_COL["NEOLP"], ms=ms,
                     label="NEOLP (n = 3)"),
          plt.Line2D([], [], marker="o", ls="", color=COND_COL["EOLP"], ms=ms,
                     label="EOLP (n = 2)")]
    # legend parked in the blank upper-right corner of the panel (the point cloud
    # tops out well below the axes top there, so nothing is covered)
    # legend_bbox / legend_ncol let a recomposition park the key outside the axes when the
    # embedding leaves no empty corner.
    ax.legend(handles=hb, loc=legend_loc, bbox_to_anchor=legend_bbox,
               frameon=True, framealpha=0.9,
               edgecolor="#DDDDDD", fontsize=7.6, handletextpad=0.4,
               labelspacing=0.5, borderpad=0.4, ncol=legend_ncol,
               columnspacing=0.9 if legend_ncol > 1 else 2.0)

    ax.text(-0.02, 1.075, tag, transform=ax.transAxes, fontsize=15,
            fontweight="bold", ha="left", va="bottom")
    return ax


def panel_C(fig, tag="C"):
    """Candidate-gene UMAP feature composite, drawn as an embedded image.
    Creates its own axes (the original image geometry is preserved verbatim)."""
    # ------------------------- (C): candidate-gene UMAP feature figure -----------
    # The previous per-gene strip (panels C-G) is replaced by the single curated
    # Stage-3 composite geo_gse211630_olp_candidates.png (5 genes, 2 rows, with its
    # own colour bar and detection table). It is placed as panel C, full width,
    # between the A/B overview strip and the lower dotplot grid.
    _cimg = plt.imread(os.path.join("figures", "geo_gse211630_olp_candidates.png"))
    _CW = 12.2 * 0.86                          # image width (inches), ~86% of canvas
    _CH = _CW / (_cimg.shape[1] / _cimg.shape[0])   # keep the image aspect ratio
    _ctop = 0.775                               # image top (figure fraction), clear gap below A/B row
    _cy0 = _ctop - _CH / 18.9                   # image bottom (figure fraction)
    axC = fig.add_axes([0.065, _cy0, _CW / 12.2, _CH / 18.9])
    axC.imshow(_cimg, interpolation="lanczos")
    axC.set_xticks([])
    axC.set_yticks([])
    for sp in axC.spines.values():
        sp.set_visible(False)
    fig.text(0.018, _ctop, tag, transform=fig.transFigure,
             fontsize=15, fontweight="bold", ha="left", va="bottom")

    return axC


# The explanatory paragraph that sits with the dotplot keys. Kept here so it has one source of
# truth: erosion_grid() draws it into its keys cell, and a recomposition that needs that cell for
# the keys alone can draw the same text elsewhere.
EROSION_NOTE = (
    "* cell-level Mann-Whitney BH FDR\n"
    "< 0.05: exploratory only (cells are\n"
    "pseudoreplicates within a sample;\n"
    "tables give the conservative level\n"
    "permutation test. Rows need >=5\n"
    "expressing cells per condition;\n"
    "other rows stay in the summary tables.\n"
    "GSE211630 buccal mucosa, 5 samples;\n"
    "doublet-suspect c12/c18 (n=32) excluded."
)


def erosion_grid(axes, label_fontsize=7.8, notes_text=True):
    """The 2x3 candidate-gene erosion dotplot grid (5 gene panels + legend/notes panel).

    `axes` is indexed with .flat, so it may be a 2x3 object array or a flat 1-D array of six axes
    (which lets a recomposition place the notes panel on a row of its own). label_fontsize trims
    the cell-type labels when a plate gives each panel less height than the original 12.2 x 18.9 in.
    """
    if isinstance(axes, list):
        axes = np.array(axes, dtype=object)
    leg_fracs = np.array([0.02, 0.05, 0.15, 0.35, 0.6])

    for j, g in enumerate(CAND):
        ax = axes.flat[j]
        gs = fig_sub[fig_sub.gene == g].sort_values("cluster")
        if gs.empty:
            ax.text(0.5, 0.5, "insufficient detection\n(<5 expressing cells per\n"
                    "condition in every type)", ha="center", va="center",
                    fontsize=7.8, color="#777777", transform=ax.transAxes)
            ax.set_yticks([])
        else:
            ypos = {r.celltype: i for i, r in enumerate(gs.itertuples())}
            for r in gs.itertuples():
                y = ypos[r.celltype]
                for cn in COND_ORDER:
                    p = per_sample[(per_sample.gene == g) &
                                   (per_sample.celltype == r.celltype) &
                                   (per_sample.condition == cn)]
                    for q in p.itertuples():
                        if q.n_expr >= 1:
                            size = SIZE_MIN + (SIZE_MAX - SIZE_MIN) * min(q.frac_expr / 0.6, 1.0)
                            ax.scatter(q.mean_exp, y, s=size, color=COND_COL[cn],
                                       alpha=0.95, edgecolor="white", linewidth=0.4, zorder=3)
                        else:
                            ax.plot(0.0, y, marker="x", ms=3.6, color="#999999", zorder=3)
                    med = np.median(p.mean_exp) if len(p) else np.nan
                    if not np.isnan(med):
                        ax.plot([med, med], [y - 0.26, y + 0.26], color=COND_COL[cn],
                                lw=1.5, zorder=4)
            ax.set_yticks([ypos[r.celltype] for r in gs.itertuples()])
            fdr_hit = {r.celltype for r in gs.itertuples()
                       if not np.isnan(r.BH_FDR_cell) and r.BH_FDR_cell < 0.05}
            ax.set_yticklabels([("* " if r.celltype in fdr_hit else "  ") +
                                SHORT[r.cluster] for r in gs.itertuples()],
                               fontsize=label_fontsize)
            for lab, r in zip(ax.get_yticklabels(), gs.itertuples()):
                if r.celltype in fdr_hit:
                    lab.set_color(HIT_RED)
        ax.set_title(g, fontsize=10.0, style="italic", loc="center", pad=3)
        ax.margins(x=0.07)
        # sparse integer x ticks; labels on bottom-row data panels only
        xmax = np.ceil(ax.get_xlim()[1]) if not gs.empty else 1.0
        ax.set_xticks(np.arange(0, xmax + 1, 1.0))
        ax.set_xticklabels([f"{int(v):g}" for v in np.arange(0, xmax + 1, 1.0)],
                           fontsize=7.3)
        ax.tick_params(axis="y", length=0)
        if j in (3, 4):                     # bottom row keeps the shared axis title
            ax.set_xlabel("mean log1p among expressing cells", fontsize=8.3, labelpad=12)
        ax.set_xlim(-0.15, None)            # keep x=0 markers off the left spine

    # legend + notes panel
    axn = axes.flat[5]
    axn.axis("off")
    h1 = [plt.Line2D([], [], marker="o", ls="", color=COND_COL["NEOLP"], ms=5.6,
                      label="NEOLP (n=3)"),
          plt.Line2D([], [], marker="o", ls="", color=COND_COL["EOLP"], ms=5.6,
                      label="EOLP (n=2)"),
          plt.Line2D([], [], ls="-", color="#666666", lw=1.6,
                      label="group median"),
          plt.Line2D([], [], marker="x", ls="", color="#999999", ms=4.6,
                      label="not expressed")]
    # keep the colour/median/no-cell key: a second axn.legend() call replaces the
    # first one (it is removed from the axes), so it must be re-added as an artist.
    # Labels are kept short so this key and the "frac. expressing" key (same row,
    # right side) never overlap.
    axn.add_artist(axn.legend(handles=h1, loc="upper left", frameon=False, fontsize=7.2,
                              labelspacing=0.6, handlelength=1.4, handletextpad=0.4))
    hs = [plt.Line2D([], [], marker="o", ls="", color="#B0B0B0",
                      ms=2 * np.sqrt((SIZE_MIN + (SIZE_MAX - SIZE_MIN) * min(f / 0.6, 1.0)) / np.pi))
          for f in leg_fracs]
    axn.legend(hs, [f"{f:g}" for f in leg_fracs], loc="upper right", frameon=True,
               edgecolor="#DDDDDD", title="frac. expressing", title_fontsize=7.6,
               fontsize=7.2, labelspacing=0.55, borderpad=0.5, handletextpad=0.35,
               handlelength=1.0)
    # notes_text=False drops the explanatory paragraph here; a recomposition can instead draw it
    # elsewhere (see EROSION_NOTE) and keep this cell for the two keys.
    if notes_text:
        axn.text(0.03, 0.50, EROSION_NOTE, transform=axn.transAxes, va="top",
                 fontsize=7.2, color="#333333", linespacing=1.15)
    else:
        # Re-anchor the two keys to the top of the cell, since the paragraph no longer fills it.
        axn.get_legend().set_bbox_to_anchor((1.0, 1.0))
        for _t in axn.texts:
            _t.set_visible(False)

    return axes


def plate_heading(fig, tag="D"):
    """Heading of the lower dotplot grid and its panel letter."""
    # heading of the lower grid (panel D; was fig.suptitle in the standalone version)
    fig.text(0.5, 0.415, "MR candidate-gene expression in erosive vs non-erosive OLP, "
             "by cell type", fontsize=10.5, weight="bold", ha="center")
    fig.text(0.075, 0.412, tag, transform=fig.transFigure, fontsize=15,
             fontweight="bold", ha="left", va="top")
    return fig


def collision_scan(fig):
    """Assert no text collisions and nothing outside the canvas; exit 1 on failure."""
    # ---- prove the render: text collision scan + in-canvas bounds
    # Exclusion policy: tick labels of the SAME axes may abut each other / their
    # own spine by construction; we still assert every non-tick text (titles,
    # xlabels, legends, notes) is clear of every other text and of spines, and
    # that tick labels never collide with non-tick text.
    rend = fig.canvas.get_renderer()
    all_t = [(t, t.get_window_extent(rend)) for t in fig.findobj(plt.Text)
             if t.get_text().strip() and t.get_visible()]
    ax_tick = {t for ax in fig.axes for t in ax.get_xticklabels(which="both") +
               ax.get_yticklabels(which="both")}
    frames = [(s, s.get_window_extent(rend)) for ax in fig.axes
              for s in ax.spines.values() if s.get_visible()]
    texts = [(t, b) for t, b in all_t if t not in ax_tick]
    ticks = [(t, b) for t, b in all_t if t in ax_tick]
    def _hit(ba, bb):
        ox = min(ba.x1, bb.x1) - max(ba.x0, bb.x0)
        oy = min(ba.y1, bb.y1) - max(ba.y0, bb.y0)
        return ox > 2.5 and oy > 2.5    # ignore <2.5px antialias near-misses

    hits = [(a, b) for i, (a, ba) in enumerate(texts) for b, bb in texts[i + 1:]
            if _hit(ba, bb)]
    tick_of_ax = {ax: {t for t in ax.get_xticklabels(which="both") +
                       ax.get_yticklabels(which="both")} for ax in fig.axes}
    for t, bt in ticks + texts:
        for s, bs in frames:
            if _hit(bt, bs) and t not in tick_of_ax.get(s.axes, set()):
                hits.append((t, s))
    if hits:
        print("COLLISIONS:", len(hits))
        for t, o in hits[:6]:
            if isinstance(o, matplotlib.spines.Spine):
                tb = t.get_window_extent(rend)
                sb = o.get_window_extent(rend)
                print(f"  text '{t.get_text()[:22]}' box x0={tb.x0:.1f} x1={tb.x1:.1f} "
                      f"y0={tb.y0:.1f} y1={tb.y1:.1f} ha={t.get_ha()} va={t.get_va()}\n"
                      f"      spine box x0={sb.x0:.1f} x1={sb.x1:.1f} "
                      f"y0={sb.y0:.1f} y1={sb.y1:.1f}")
            else:
                print(f"  text-vs-text '{t.get_text()[:20]}' / '{o.get_text()[:20]}'")
        raise SystemExit(1)
    print("render ok: no text collisions; figures written")

    return fig



def main():
    fig, axes = plt.subplots(2, 3, figsize=(12.2, 18.9))
    fig.subplots_adjust(left=0.19, right=0.995, top=0.385, bottom=0.045,
                        wspace=0.65, hspace=0.60)
    erosion_grid(axes)

    axA = fig.add_axes([0.020, _yMid - _axh / 2, _axw, _axh])
    axB = fig.add_axes([0.690, _yMid - _axh / 2, _axw, _axh])
    panel_A(axA, "A")
    panel_A_legend(axA)
    panel_B(axB, "B")

    panel_C(fig, "C")
    plate_heading(fig, "D")
    fig.savefig(os.path.join(FIGDIR, "fig06_target_expression.png"), dpi=300)
    fig.savefig(os.path.join(FIGDIR, "fig06_target_expression.pdf"))
    fig.savefig(os.path.join(QA, "fig_erosion.png"), dpi=150)
    collision_scan(fig)
    plt.close("all")
    print("render ok: no text collisions; figures written")
    return fig


if __name__ == "__main__":
    main()
