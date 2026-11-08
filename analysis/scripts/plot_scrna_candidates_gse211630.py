# -*- coding: utf-8 -*-
"""Stage 3 figures & cell-type annotation aid for GSE211630 scanpy objects.

Outputs (figures/):
  geo_gse211630_olp_umap.png         OLP-only: Leiden UMAP | condition UMAP
  geo_gse211630_olp_candidates.png   OLP-only: 5 candidate genes on UMAP
  geo_gse211630_olp_marker_dotplot.png  OLP-only: cluster x lineage-marker dotplot
  geo_gse211630_all_umap_samples.png  all-samples: sample-facet UMAP + 3 candidates
Plus results/tables/geo_gse211630_candidate_loc.csv and a celltype suggestion
JSON (data/processed/geo_GSE211630/scanpy/celltype_suggestion.json).
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

# The figure-style skill shipped `kernel.py` when this script was written; it now ships
# `runtime.py` (renamed 2026-09-21) with the same apply_figure_style/set_frame helpers. Resolve
# whichever is present so the script keeps running after the rename.
_kernel = r"\\?\D:\wisp-science\skills\figure-style\kernel.py"
_runtime = r"\\?\D:\wisp-science\skills\figure-style\runtime.py"
sidecar = _kernel if os.path.exists(_kernel) else _runtime
exec(compile(open(sidecar, encoding="utf-8").read(), "kernel.py", "exec"))
apply_figure_style(sizes=(8, 7, 6))

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
os.chdir(ROOT)
OUT = os.path.join("data", "processed", "geo_GSE211630")
SCAN = os.path.join(OUT, "scanpy")
FIGDIR = "figures"
os.makedirs(FIGDIR, exist_ok=True)

CAND = ["GDF11", "SUOX", "EGFL8", "TAPBP", "HLA-DPA1"]
OLP_COL = {"NEOLP": "#3B6FB5", "EOLP": "#C26401"}

print("loading analysed objects...", flush=True)
olp = sc.read_h5ad(os.path.join(SCAN, "gse211630_olp_lognorm.h5ad"))
allobj = sc.read_h5ad(os.path.join(SCAN, "gse211630_all_lognorm.h5ad"))
leid = "leiden_0.5"
print(f"olp {olp.shape} clusters={olp.obs[leid].nunique()} | all {allobj.shape} clusters={allobj.obs[leid].nunique()}", flush=True)

umap = np.asarray(olp.obsm["X_umap"])
umap_all = np.asarray(allobj.obsm["X_umap"])
rng = np.random.default_rng(0)


def umap_ax(ax, coords, col, cmap=None, vmin=None, vmax=None,
            categorical=False, cmap_cat=None, title="", s=0.35, alpha=0.6,
            legend_ax=None):
    ax.scatter(coords[:, 0], coords[:, 1], c=col, cmap=cmap, vmin=vmin, vmax=vmax,
               s=s, alpha=alpha, linewidths=0, rasterized=True)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(title, loc="left", fontsize=8, weight="bold")


# ---------------- fig 1: OLP UMAP (Leiden + condition) ----------------
# Right margin is reserved for both legends (axes-anchored legends overflowed
# the canvas and clipped cluster labels 11-20).
fig, axes = plt.subplots(1, 2, figsize=(12.6, 4.4))
fig.subplots_adjust(left=0.01, right=0.79, top=0.93, bottom=0.02, wspace=0.10)
clusters = sorted(olp.obs[leid].astype(int).unique())
cols_cl = pd.Series(olp.obs[leid]).astype(int).map(
    {c: plt.cm.tab20(c % 20) for c in clusters})
umap_ax(axes[0], umap, cols_cl, title=f"OLP Leiden (res 0.5), n={olp.n_obs:,}")
cl_hand = [plt.Line2D([], [], marker="o", ls="", color=plt.cm.tab20(c % 20),
                       label=str(c), alpha=0.85) for c in clusters]
fig.legend(handles=cl_hand, loc="center left", bbox_to_anchor=(0.815, 0.60),
           ncol=2, frameon=False, fontsize=6, handletextpad=0.4,
           columnspacing=0.8, borderpad=0.2)
umap_ax(axes[1], umap, olp.obs["condition"].map(OLP_COL),
        title="Condition (NEOLP vs EOLP)")
hand = [plt.Line2D([], [], marker="o", ls="", color=c, label=k, alpha=0.85)
        for k, c in OLP_COL.items()]
fig.legend(handles=hand, loc="center left", bbox_to_anchor=(0.815, 0.14),
           ncol=2, frameon=False, fontsize=7, handletextpad=0.4)
fig.savefig(os.path.join(FIGDIR, "geo_gse211630_olp_umap.png"), dpi=300)
plt.close(fig)
print("saved geo_gse211630_olp_umap.png", flush=True)

# ---------------- fig 2: candidate genes on OLP UMAP ----------------
# Shared colour scale across all panels (0-3 log1p) so panels are directly
# comparable; cells with zero expression are drawn light-grey ("not detected")
# instead of near-black, and expressing cells are overlaid in colour at larger
# size so sparse genes (GDF11/SUOX/EGFL8) remain visible.
def panel_candidate_umap(ax, g):
    """One candidate gene on the OLP UMAP: grey landscape, expressing cells in magma.

    Returns (frac, max, n_expr, per_condition) or None when the gene is absent from the object.
    Shared 0-VMAX_EXPR scale and grey no-detection backdrop keep the panels comparable, and
    expressing cells are overlaid larger so sparse genes stay visible.
    """
    if g not in olp.var_names:
        return None
    v = np.asarray(olp[:, g].X.todense()).ravel() if hasattr(olp[:, g].X, "todense") \
        else np.asarray(olp[:, g].X).ravel()
    # grey landscape: every cell (no-detection cells invisible against it)
    ax.scatter(umap[:, 0], umap[:, 1], c=GRAY_NODET, s=0.28, alpha=0.45,
               linewidths=0, rasterized=True)
    m = v > 0
    if m.any():
        ax.scatter(umap[m, 0], umap[m, 1], c=v[m], cmap="magma",
                   vmin=0, vmax=VMAX_EXPR, s=0.85, alpha=0.85,
                   linewidths=0, rasterized=True)
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(g, loc="left", fontsize=8, weight="bold")
    per = {}
    for cond in ("NEOLP", "EOLP"):
        mm = olp.obs["condition"] == cond
        per[cond] = {"frac": float((v[mm] > 0).mean()),
                     "mean_log": float(v[mm].mean())}
    return float((v > 0).mean()), float(v.max()), int((v > 0).sum()), per


ncol = 3
nrow = 2
VMAX_EXPR = 3.0  # log1p(norm. counts/10k); values >3 saturate (rare tails)
GRAY_NODET = "#C7C7C7"
fig, axes = plt.subplots(nrow, ncol, figsize=(11.2, 6.8))
fig.subplots_adjust(left=0.01, right=0.90, top=0.95, bottom=0.01,
                    wspace=0.06, hspace=0.20)
expr_frac, summary_rows, expr_max = {}, [], {}
for i, g in enumerate(CAND):
    ax = axes.flat[i]
    res = panel_candidate_umap(ax, g)
    if res is not None:
        frac, vmax_, nexp, per = res
        expr_frac[g] = frac
        expr_max[g] = vmax_
        summary_rows.append({"gene": g, "present_in_olp": True,
                             "n_cells_expressing": nexp,
                             "frac_cells": round(frac, 4),
                             **{f"{c}_frac": round(per[c]["frac"], 4) for c in per},
                             **{f"{c}_mean_log1p": round(per[c]["mean_log"], 4) for c in per}})
    else:
        ax.text(0.5, 0.5, f"{g}\n(row absent in OLP object)", ha="center",
                va="center", fontsize=7, color="#888888")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(g, loc="left", fontsize=8, weight="bold")
        summary_rows.append({"gene": g, "present_in_olp": False,
                             "n_cells_expressing": None, "frac_cells": None})
# 6th cell: shared-scale note + per-gene stats table (was blank);
# shared colorbar lives in the reserved right margin (no tick/text collisions).
axn = axes.flat[len(CAND)]
axn.axis("off")
_note_lines = ("grey = no detection\n"
               "scale shared 0-3 log1p(norm. counts/10k);\n"
               "values >3 saturate (rare tails)\n\n"
               "OLP, n = 46,478 cells\n\n"
               "gene      %expr   max\n" +
               "\n".join(
                   f"{g:<9} {expr_frac[g]*100:5.1f}%   {expr_max[g]:4.1f}"
                   for g in CAND if g in olp.var_names))
axn.text(0.02, 0.985, _note_lines, ha="left", va="top", fontsize=6.4,
         color="#222222", linespacing=1.55, family="monospace")
_cax = fig.add_axes([0.925, 0.13, 0.022, 0.74])
cb = fig.colorbar(matplotlib.cm.ScalarMappable(
    norm=matplotlib.colors.Normalize(vmin=0, vmax=VMAX_EXPR), cmap="magma"),
    cax=_cax)
cb.set_label("log1p(norm. counts/10k)", fontsize=6.5)
cb.ax.tick_params(labelsize=5.5)
fig.savefig(os.path.join(FIGDIR, "geo_gse211630_olp_candidates.png"), dpi=300)
plt.close(fig)
print("saved geo_gse211630_olp_candidates.png", flush=True)
print("expr fraction:", expr_frac, flush=True)

# ---------------- fig 3: marker dotplot + celltype suggestion ----------------
MARKERS = {
    "T cell": ["CD3D", "CD3E", "TRAC", "IL7R", "CD2"],
    "NK": ["NKG7", "GNLY", "KLRD1", "KLRB1"],
    "B cell": ["MS4A1", "CD79A", "CD79B", "CD19"],
    "Plasma cell": ["MZB1", "SDC1", "JCHAIN", "DERL3", "XBP1"],
    "Myeloid": ["LYZ", "C1QA", "C1QB", "FCGR3A", "CD68", "HLA-DRA"],
    "Mast": ["TPSAB1", "TPSB2", "CPA3", "MS4A2"],
    "Epith basal": ["KRT5", "KRT14", "KRT17", "TP63", "COL17A1"],
    "Epith diff": ["KRT13", "KRT4", "KRT6A", "KRT16"],
    "Fibroblast": ["COL1A1", "COL3A1", "DCN", "LUM", "COL6A1"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5"],
    "Pericyte/SMC": ["ACTA2", "TAGLN", "MYH11", "RGS5"],
    "Proliferating": ["MKI67", "TOP2A", "STMN1"],
}
present_m = {grp: [g for g in genes if g in olp.var_names] for grp, genes in MARKERS.items()}
absent = {grp: [g for g in genes if g not in olp.var_names] for grp, genes in MARKERS.items()}
print("marker genes absent from OLP object:", {k: v for k, v in absent.items() if v}, flush=True)

X = olp[:, sum(present_m.values(), [])].X  # sparse log1p
X = X.tocsc() if hasattr(X, "tocsc") else X
lab = olp.obs[leid].astype(int).values
clu_list = sorted(np.unique(lab))
gene_order, gene_grp = [], []
mean_exp = np.zeros((len(clu_list), X.shape[1]))
frac_exp = np.zeros((len(clu_list), X.shape[1]))
for ci, c in enumerate(clu_list):
    sub = X[lab == c]
    mean_exp[ci] = np.asarray(sub.mean(axis=0)).ravel()
    frac_exp[ci] = np.asarray((sub > 0).mean(axis=0)).ravel()
for grp, genes in present_m.items():
    for g in genes:
        gene_order.append(g)
        gene_grp.append(grp)
gi = {g: k for k, g in enumerate(sum(present_m.values(), []))}
z = mean_exp.copy()
for j in range(z.shape[1]):
    sd = z[:, j].std()
    z[:, j] = (z[:, j] - z[:, j].mean()) / (sd if sd > 0 else 1.0)

fig, ax = plt.subplots(figsize=(14.6, 4.7))
fig.subplots_adjust(left=0.19, right=0.86, top=0.93, bottom=0.12)
grp_bounds = []
for grp, genes in present_m.items():
    idx = [gi[g] for g in genes]
    grp_bounds.append((grp, min(idx), max(idx)))
pos = np.arange(len(gene_order))
# grid lines separating groups
for grp, lo, hi in grp_bounds:
    ax.axvline(lo - 0.5, color="#DDDDDD", lw=0.5)
scat = ax.scatter(np.repeat(pos[None, :], len(clu_list), axis=0),
                  np.repeat(np.array(clu_list)[:, None], len(gene_order), axis=1),
                  c=z, cmap="RdBu_r", vmin=-max(1.0, np.abs(z).max() * 0.8) if z.size else -1,
                  vmax=max(1.0, np.abs(z).max() * 0.8) if z.size else 1,
                  s=(2.0 + 78.0 * np.clip(frac_exp, 0, 1) ** 0.5).ravel(),
                  linewidths=0, alpha=0.9, rasterized=True)
ax.set_xticks(pos)
ax.set_xticklabels(gene_order, rotation=90, fontsize=5.2)
ax.set_yticks(clu_list)
ax.set_ylabel("Leiden cluster", fontsize=7)
ax.set_ylim(clu_list[0] - 0.5, clu_list[-1] + 0.5)
ax.set_xlim(-0.5, len(gene_order) - 0.5)
# top group labels
for grp, lo, hi in grp_bounds:
    ax.text((lo + hi) / 2, len(clu_list) - 0.35, grp, ha="center", va="bottom",
            fontsize=6, color="#444444")
cbar = fig.colorbar(scat, ax=ax, fraction=0.02, pad=0.01)
cbar.set_label("z(mean log1p expr.)", fontsize=6)
cbar.ax.tick_params(labelsize=5)
# dot-size legend (fraction of cells expressing each marker in each cluster)
_size_prox = []
for _f in (0.05, 0.25, 0.5, 0.75, 1.0):
    _s = 2.0 + 78.0 * np.sqrt(_f)
    _size_prox.append(plt.Line2D([], [], marker="o", ls="", color="#8C8C8C",
                                 markersize=2 * np.sqrt(_s / np.pi),
                                 label=f"{_f:g}"))
fig.legend(handles=_size_prox, loc="center left", bbox_to_anchor=(0.895, 0.5),
           frameon=False, fontsize=7, title="frac. expressing",
           title_fontsize=7.5, labelspacing=0.65, borderpad=0.2, handletextpad=0.35)
fig.savefig(os.path.join(FIGDIR, "geo_gse211630_olp_marker_dotplot.png"), dpi=300)
plt.close(fig)
print("saved geo_gse211630_olp_marker_dotplot.png", flush=True)

# cell-type suggestion: argmax group z per cluster with margin
score = {}
for grp in present_m:
    if not present_m[grp]:
        continue
    score[grp] = z[:, [gi[g] for g in present_m[grp]]].mean(axis=1)
S = pd.DataFrame(score, index=clu_list)
sugg = {}
for c in clu_list:
    top = S.loc[c].sort_values(ascending=False)
    sugg[str(c)] = {"top1": top.index[0], "score_top1": round(float(top.iloc[0]), 2),
                    "top2": top.index[1], "score_top2": round(float(top.iloc[1]), 2),
                    "margin": round(float(top.iloc[0] - top.iloc[1]), 2)}
with open(os.path.join(SCAN, "celltype_suggestion.json"), "w") as f:
    json.dump({"leiden_col": leid, "per_cluster": sugg, "markers_absent": absent},
              f, indent=2)
print("celltype suggestion (top1 | top2 | margin):", flush=True)
for c in clu_list:
    s = sugg[str(c)]
    print(f"  cluster {c}: {s['top1']} ({s['score_top1']}) vs {s['top2']} ({s['score_top2']}), margin {s['margin']}", flush=True)

# ---------------- fig 4a: all-samples sample-facet UMAP ----------------
# Reading order normal -> NEOLP -> EOLP (severity progression; normal as control).
_SAMPLE_ORDER = ["normal-1", "NEOLP-1", "NEOLP-2", "NEOLP-3", "EOLP-1", "EOLP-2"]
_present = list(allobj.obs["sample"].unique())
assert set(_SAMPLE_ORDER) == set(_present), (_SAMPLE_ORDER, _present)
samples_all = _SAMPLE_ORDER
fig, axes = plt.subplots(2, 3, figsize=(10.8, 7.0))
fig.subplots_adjust(left=0.01, right=0.99, top=0.96, bottom=0.01,
                    wspace=0.05, hspace=0.24)
for si, samp in enumerate(samples_all):
    ax = axes.flat[si]
    m = allobj.obs["sample"] == samp
    sub = umap_all[m]
    ax.scatter(sub[:, 0], sub[:, 1], c="#7A7A7A", s=0.3, alpha=0.45,
               linewidths=0, rasterized=True)
    ax.set_title(f"{samp} (n={int(m.sum()):,})", loc="left", fontsize=7.5,
                 weight="bold")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
fig.savefig(os.path.join(FIGDIR, "geo_gse211630_all_umap_samples.png"), dpi=300)
plt.close(fig)
print("saved geo_gse211630_all_umap_samples.png", flush=True)

# ---------------- fig 4b: candidates on all-samples object ----------------
cand_all = [g for g in CAND if g in allobj.var_names]  # EGFL8/TAPBP/HLA-DPA1
miss_all = [g for g in CAND if g not in allobj.var_names]  # GDF11/SUOX
fig, axes = plt.subplots(2, 2, figsize=(8.6, 8.6))
fig.subplots_adjust(left=0.01, right=0.88, top=0.95, bottom=0.01,
                    wspace=0.08, hspace=0.20)
expr_frac_all = {}
for i, g in enumerate(cand_all):
    ax = axes.flat[i]
    v = np.asarray(allobj[:, g].X.todense()).ravel() if hasattr(allobj[:, g].X, "todense") \
        else np.asarray(allobj[:, g].X).ravel()
    ax.scatter(umap_all[:, 0], umap_all[:, 1], c=GRAY_NODET, s=0.25, alpha=0.45,
               linewidths=0, rasterized=True)
    m = v > 0
    if m.any():
        ax.scatter(umap_all[m, 0], umap_all[m, 1], c=v[m], cmap="magma",
                   vmin=0, vmax=VMAX_EXPR, s=0.8, alpha=0.85,
                   linewidths=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    frac_all = (v > 0).mean()
    expr_frac_all[g] = float(frac_all)
    ax.set_title(f"{g} (all 6 samples)", loc="left", fontsize=7.5, weight="bold")
ax = axes.flat[len(cand_all)]
ax.axis("off")
_note = (f"grey = no detection; scale shared 0-{VMAX_EXPR:g} log1p "
         f"(norm. counts/10k)\n\n"
         "Row absent (all-object gene universe =\nnormal-1 12,021 genes; "
         "Stage-1 check):\n" + "\n".join(f"• {g}" for g in miss_all) + "\n\n" +
         "EGFL8: row present but sparse — " +
         "expressing cells " + f"{expr_frac_all.get('EGFL8', 0)*100:.1f}%")
ax.text(0.04, 0.97, _note, ha="left", va="top", fontsize=7.2,
        color="#666666", transform=ax.transAxes, linespacing=1.5)
_cax = fig.add_axes([0.905, 0.12, 0.022, 0.76])
_cb2 = fig.colorbar(matplotlib.cm.ScalarMappable(
    norm=matplotlib.colors.Normalize(vmin=0, vmax=VMAX_EXPR), cmap="magma"),
    cax=_cax)
_cb2.set_label("log1p(norm. counts/10k)", fontsize=6.5)
_cb2.ax.tick_params(labelsize=5.5)
fig.savefig(os.path.join(FIGDIR, "geo_gse211630_all_candidates.png"), dpi=300)
plt.close(fig)
print("saved geo_gse211630_all_candidates.png", flush=True)

# ---------------- candidate summary CSV ----------------
dfc = pd.DataFrame(summary_rows)
# enrich with raw-max-UMI from Stage-1 check
with open(os.path.join(OUT, "candidate_genes_check.json"), encoding="utf-8") as f:
    cgc = json.load(f)
olp_gsm = ["GSM6481638", "GSM6481639", "GSM6481640", "GSM6481641", "GSM6481642"]
rawmax = {}
for g in CAND:
    vals = [cgc[f"{gs}_count_data-{nm}.txt.gz"]["genes"].get(g, {}).get("max_umi", 0)
            for gs, nm in zip(olp_gsm, ["NEOLP-1", "NEOLP-2", "NEOLP-3", "EOLP-1", "EOLP-2"])]
    rawmax[g] = max(vals) if vals else 0
dfc["raw_max_umi_olp"] = dfc["gene"].map(rawmax)
dfc.to_csv(os.path.join("results", "tables", "geo_gse211630_candidate_loc.csv"),
           index=False)
print(dfc.to_string(index=False), flush=True)
print("DONE", flush=True)
