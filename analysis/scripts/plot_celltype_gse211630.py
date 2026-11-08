# -*- coding: utf-8 -*-
"""Stage 4 (option B): publication figures for formal cell-type annotation.
Fig 1: annotated UMAP (fine cell types, hue-family palette; doublet-suspect grey).
Fig 2: candidate genes x annotated cell type dotplot (size=frac expressing;
       colour=mean log1p among expressing; ring = coarse-family Fisher FDR<0.05).
Outputs -> figures/geo_gse211630_olp_celltype_umap.png,
            figures/geo_gse211630_candidate_celltype_dotplot.png
"""
import colorsys
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
os.chdir(ROOT)
SCAN = os.path.join("data", "processed", "geo_GSE211630", "scanpy")
TAB = os.path.join("results", "tables")
FIGDIR = "figures"
QA = os.path.join("analysis", "scrna", "_qa")
os.makedirs(QA, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7,
                     "axes.linewidth": 0.6, "svg.fonttype": "none"})

# ------------------------------------------------------------------ annotation
olp = sc.read_h5ad(os.path.join(SCAN, "gse211630_olp_lognorm.h5ad"))
olp.obs["leiden_0.5"] = olp.obs["leiden_0.5"].astype(int)
ob = pd.read_csv(os.path.join(SCAN, "gse211630_olp_obs_annotation.csv"), index_col=0)
assert (olp.obs.index == ob.index).all()
olp.obs["celltype_sub"] = pd.Categorical(ob["celltype_sub"].values)
olp.obs["celltype_family"] = pd.Categorical(ob["celltype_family"].values)
olp.obs["doublet_suspect"] = ob["doublet_suspect"].values

FAMILY_ORDER = ["T/NK", "B/Plasma", "Myeloid", "Mast", "Epithelial",
                "Fibroblast", "Pericyte/SMC", "Endothelial"]
CL_SUB = {  # cluster -> sub annotation (mirrors scrna_annotate_gse211630.py)
    0: "Fibroblast (SFRP2+/COL1A2+)", 1: "Fibroblast (CXCL14+/APOD+)",
    2: "Pericyte/SMC",
    3: "T cell (resting/recirculating, CD4/CD8)",
    4: "T/NK (stressed, HSP-high)", 5: "T/NK-like (low-complexity)",
    6: "CD8+ T (GZMK+ effector-memory)", 7: "CD4+ Treg (FOXP3+)",
    8: "T cell (naive-like, IL7Rhigh)",
    9: "Plasmablast (GZMB+)", 10: "Macrophage (LYZ+/HLA-DR+)",
    11: "Epithelial keratinocyte (basal/proliferating)",
    13: "Mast cell", 14: "B cell (MS4A1+)", 15: "Plasma cell (SDC1+)",
    16: "CD8+ T (GZMB+/CXCL13+ cytotoxic-exhausted-like)",
    17: "Endothelial (lymphatic, CCL21+)", 19: "Endothelial (blood, ENG+/ADGRL4+)",
    20: "Myeloid (S100-high inflammatory)"}
CL_FAM = {0: "Fibroblast", 1: "Fibroblast", 2: "Pericyte/SMC", 3: "T/NK", 4: "T/NK",
          5: "T/NK", 6: "T/NK", 7: "T/NK", 8: "T/NK", 9: "B/Plasma", 10: "Myeloid",
          11: "Epithelial", 13: "Mast", 14: "B/Plasma", 15: "B/Plasma", 16: "T/NK",
          17: "Endothelial", 19: "Endothelial", 20: "Myeloid"}
DOUBLET = 12, 18

# per-family HLS palette (hue coherence within family; colourblind-tolerant)
FAM_HLS = {"T/NK": (0.00, 0.47, 0.92), "B/Plasma": (0.585, 0.44, 0.80),
           "Myeloid": (0.34, 0.47, 0.72), "Mast": (0.80, 0.52, 0.62),
           "Epithelial": (0.105, 0.42, 0.88), "Fibroblast": (0.145, 0.40, 0.72),
           "Pericyte/SMC": (0.07, 0.26, 0.55), "Endothelial": (0.50, 0.42, 0.72)}
# 7 T/NK subtypes need a full HLS override: pure-red lightness rungs saturate
# (red=max at L>~0.52 for high S), so light members must desaturate to read as
# distinct pinks/salmons while dark members stay saturated (crimson -> orange-red)
_SUB_HLS = {
    3: (0.000, 0.50, 0.95),   # pure red
    4: (0.000, 0.66, 0.50),   # coral-pink (clearly lighter/pinker than c3)
    5: (0.015, 0.82, 0.28),   # pale pink (lightest)
    6: (0.985, 0.24, 0.80),   # burgundy (darkest)
    7: (0.035, 0.55, 0.90),   # orange-red
    8: (0.050, 0.74, 0.42),   # salmon (vs c5 pale pink)
    16: (0.000, 0.36, 0.88),  # dark red
}
_SUB_LADJ = {9: +0.24, 14: 0.00, 15: -0.15,
             10: 0.00, 20: +0.15, 11: 0.00, 0: 0.00, 1: +0.20, 2: 0.00,
             17: +0.25, 19: 0.00}


def sub_color(cl):
    if cl in _SUB_HLS:
        h, l2, s = _SUB_HLS[cl]
        return colorsys.hls_to_rgb(h, l2, s)
    h, l, s = FAM_HLS[CL_FAM[cl]]
    l2 = float(np.clip(l + _SUB_LADJ.get(cl, 0.0), 0.16, 0.82))
    return colorsys.hls_to_rgb(h, l2, s)


sub_legend_order = []
for fam in FAMILY_ORDER:
    for cl in sorted(CL_FAM):
        if CL_FAM[cl] == fam:
            sub_legend_order.append((cl, CL_SUB[cl], fam))
colmap = {cl: sub_color(cl) for cl in CL_FAM}
colmap.update({12: (0.76, 0.76, 0.76), 18: (0.76, 0.76, 0.76)})

umap = np.asarray(olp.obsm["X_umap"])
leid = olp.obs["leiden_0.5"].values.astype(int)

# ====================================================================== FIG 1
fig, ax = plt.subplots(figsize=(9.2, 7.2))
fig.subplots_adjust(left=0.005, right=0.66, top=0.97, bottom=0.005)
n_tot = len(olp)
# base canvas first: doublet-suspect drawn UNDER, then re-drawn on top (grey,
# visible markers) since they are excluded from stats but must remain locatable
for dcl in DOUBLET:
    m = leid == dcl
    ax.scatter(umap[m, 0], umap[m, 1], s=1.1, c=[colmap[dcl]] * int(m.sum()),
               linewidths=0, alpha=0.8, rasterized=True)
for cl in sorted(CL_FAM):
    m = leid == cl
    ax.scatter(umap[m, 0], umap[m, 1], s=1.4, c=[colmap[cl]] * int(m.sum()),
               linewidths=0, alpha=0.82, rasterized=True)
for dcl in DOUBLET:
    m = leid == dcl
    ax.scatter(umap[m, 0], umap[m, 1], s=3.6, c=[colmap[dcl]] * int(m.sum()),
               linewidths=0.4, edgecolors="#444444", alpha=0.98, rasterized=True)
ax.set_xticks([])
ax.set_yticks([])
for sp in ax.spines.values():
    sp.set_visible(False)
ax.set_title("OLP (5 samples, n=46,478): 21 Leiden clusters \u2192 19 annotated "
             "cell types (8 families)", loc="left", fontsize=8, weight="bold")

handles, labels = [], []
for fam in FAMILY_ORDER:
    handles.append(plt.Line2D([], [], ls="-", color="#666666", lw=1.4, label=""))
    labels.append(fam)
    for cl, sub, fam2 in sub_legend_order:
        if fam2 != fam:
            continue
        handles.append(plt.Line2D([], [], marker="o", ls="", markersize=6.8,
                                  color=colmap[cl], alpha=0.95))
        labels.append(f"c{cl}  {sub}")
handles.append(plt.Line2D([], [], marker="o", ls="", markersize=5.6,
                          color=(0.76, 0.76, 0.76), markeredgecolor="#444444",
                          markeredgewidth=0.6))
labels.append("doublet-suspect (c12/c18, n=32, excluded from stats)")
lg = fig.legend(handles, labels, loc="center left", bbox_to_anchor=(0.672, 0.5),
                frameon=True, framealpha=0.95, edgecolor="#DDDDDD", fontsize=6.4,
                labelspacing=0.42, handletextpad=0.45, borderpad=0.5)
# bold family headers
for t, fam in zip(lg.get_texts(), labels):
    if fam in FAMILY_ORDER:
        t.set_weight("bold")
        t.set_fontsize(6.6)
fig.savefig(os.path.join(FIGDIR, "geo_gse211630_olp_celltype_umap.png"), dpi=300)
fig.savefig(os.path.join(QA, "fig1.png"), dpi=150)

# ====================================================================== FIG 2
cand = ["GDF11", "SUOX", "EGFL8", "TAPBP", "HLA-DPA1"]
ann = pd.read_csv(os.path.join(TAB, "geo_gse211630_celltype_annotation.csv"))
sub_of = dict(zip(ann["leiden_0.5"], ann["celltype_sub"]))
# short row labels (c# = Leiden cluster; cross-ref with annotation table)
SHORT = {0: "c0 Fibroblast SFRP2+", 1: "c1 Fibroblast CXCL14+", 2: "c2 Pericyte/SMC",
         3: "c3 T resting/recirculating", 4: "c4 T/NK stressed (HSP+)",
         5: "c5 T/NK low-complexity", 6: "c6 CD8+ T GZMK+ EM",
         7: "c7 CD4+ Treg", 8: "c8 T naive-like (IL7R+)",
         9: "c9 Plasmablast (GZMB+)", 10: "c10 Macrophage (LYZ+)",
         11: "c11 Epithelial basal/prolif.",
         13: "c13 Mast", 14: "c14 B cell", 15: "c15 Plasma cell",
         16: "c16 CD8+ T GZMB+/CXCL13+", 17: "c17 Endo lymphatic (CCL21+)",
         19: "c19 Endo blood (ENG+)", 20: "c20 Myeloid S100-hi"}
cl_list = sorted(CL_FAM)
Xc = np.asarray(olp[:, cand].X.todense()) if hasattr(olp[:, cand].X, "todense") \
    else np.asarray(olp[:, cand].X)

rows = []
for cl in cl_list:
    m = leid == cl
    for j, g in enumerate(cand):
        v = Xc[m, j]
        ne = int((v > 0).sum())
        rows.append({"cluster": cl, "sub": sub_of[cl], "family": CL_FAM[cl],
                     "gene": g, "frac": ne / m.sum(),
                     "mean_expr": float(v[v > 0].mean()) if ne else 0.0})
dp = pd.DataFrame(rows)

# coarse-family Fisher FDR for the ring annotation
enr = pd.read_csv(os.path.join(TAB, "geo_gse211630_candidate_celltype_enrichment.csv"))
sig = set(zip(enr["gene"][enr["fisher_FDR_BH"] < 0.05],
              enr["family"][enr["fisher_FDR_BH"] < 0.05]))

fig, ax = plt.subplots(figsize=(13.8, 7.8))
fig.subplots_adjust(left=0.175, right=0.815, top=0.95, bottom=0.085)
ypos = {cl: i for i, cl in enumerate(reversed(cl_list))}
# family separators
rows_by_y = sorted([(ypos[cl], CL_FAM[cl], cl) for cl in cl_list])
prev = None
for y, fam, cl in rows_by_y:
    if prev is not None and fam != prev:
        ax.axhline(y - 0.5, color="#BBBBBB", lw=0.9)
    prev = fam

VMAX = 4.0
_fam2dark = {f: colorsys.hls_to_rgb(FAM_HLS[f][0],
                                    min(0.34, max(0.23, FAM_HLS[f][1] * 0.72)),
                                    max(0.55, FAM_HLS[f][2])) for f in FAMILY_ORDER}
for r in dp.itertuples():
    x = cand.index(r.gene)
    y = ypos[r.cluster]
    ring = (r.gene, r.family) in sig
    area = 8.0 + 620.0 * float(r.frac) ** 0.8          # marker area (pt^2)
    ax.scatter(x, y, s=area,
               c=[r.mean_expr], cmap="magma", vmin=0, vmax=VMAX,
               linewidths=1.0 if ring else 0, edgecolors="#111111" if ring else "none",
               alpha=0.94)
ax.set_xticks(range(len(cand)))
ax.set_xticklabels(cand, fontsize=8.5, style="italic")
ax.set_yticks([ypos[cl] for cl in cl_list])
ax.set_yticklabels([SHORT[cl] for cl in cl_list], fontsize=7)
for t in ax.get_yticklabels():
    cl = cl_list[[ypos[c] for c in cl_list].index(t.get_y() - 0 if False else ypos[[c for c in cl_list if SHORT[c] == t.get_text()][0]])]
    t.set_color(_fam2dark[CL_FAM[cl]])
ax.tick_params(axis="y", length=0)
ax.set_ylim(-0.6, len(cl_list) - 0.4)
ax.set_xlim(-0.6, len(cand) - 0.5)
ax.set_title("Candidate gene expression across annotated OLP cell types "
             "(n=46,446; 19 types)", loc="left", fontsize=8.5, weight="bold")
ax.set_ylabel("cell type (Leiden cluster)", fontsize=7.5)

cbar = fig.colorbar(plt.cm.ScalarMappable(
    norm=matplotlib.colors.Normalize(vmin=0, vmax=VMAX), cmap="magma"),
    ax=ax, fraction=0.018, pad=0.012)
cbar.set_label("mean log1p expr. (expressing cells)", fontsize=7)
cbar.ax.tick_params(labelsize=6)
# size legend (frac expressing) - markersize = dot diameter in points
_proxy = []
for f in (0.02, 0.1, 0.25, 0.5, 0.85):
    area = 8.0 + 620.0 * f ** 0.8
    _proxy.append(plt.Line2D([], [], marker="o", ls="", color="#B0B0B0",
                             markersize=2 * np.sqrt(area / np.pi), label=f"{f:g}"))
# largest proxy circle is d~26.5pt -> keep row pitch well above it
# so circles do not touch (labelspacing is in units of fontsize)
fig.legend(handles=_proxy, loc="center left", bbox_to_anchor=(0.845, 0.30),
           frameon=True, framealpha=0.95, edgecolor="#DDDDDD", fontsize=8,
           title="frac. expressing", title_fontsize=8, labelspacing=3.4,
           borderpad=0.5, handletextpad=0.4)
_ring = plt.Line2D([], [], marker="o", ls="", markersize=7, color="#DDDDDD",
                   markeredgecolor="#111111", markeredgewidth=0.8)
fig.legend(handles=[_ring],
           labels=["ring = family-level Fisher FDR < 0.05"],
           loc="center left", bbox_to_anchor=(0.845, 0.90), frameon=False,
           fontsize=6.6, handletextpad=0.4)
_fam_txt = "  ".join(f"{fam}" for fam in FAMILY_ORDER)
fig.text(0.175, 0.012, "row label colour = cell-type family (red T/NK | blue "
         "B/Plasma | green Myeloid | purple Mast | amber Epithelial | olive "
         "Fibroblast | brown Pericyte/SMC | teal Endothelial); "
         "doublet-suspect clusters c12/c18 (n=32) excluded",
         ha="left", va="bottom", fontsize=6.2, color="#555555")
fig.savefig(os.path.join(FIGDIR, "geo_gse211630_candidate_celltype_dotplot.png"),
            dpi=300)
fig.savefig(os.path.join(QA, "fig2.png"), dpi=150)
print("figures written", flush=True)
