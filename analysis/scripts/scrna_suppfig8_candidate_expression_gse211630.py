# -*- coding: utf-8 -*-
"""Supplementary Figure 8 - where the five candidate genes are expressed in OLP.

New plate added when the supplementary set was reordered to follow the
manuscript narrative (2026-09-13): the expression landscape of the candidate
genes now follows the marker panel (Supplementary Figure 7) and precedes the
proportion and cell-cell communication figures.

Restructured 2026-09-30 for the main-figure re-plan: each panel is now a function that draws into
a caller-supplied rectangle, so panels A/B/C can be recomposed into new main figures while this
script's own output stays unchanged. `main()` still passes the original rectangles, so it
reproduces supp_fig08_candidate_expression.{png,pdf} exactly.

Panels
  A  candidate gene x annotated cell type dot plot. Dot colour is the mean
     log1p expression among expressing cells; dot AREA is proportional to the
     fraction of cells expressing the gene (values at or below 0.01 are drawn at
     a small floor size so that a sparse detection stays visible, and a cell type
     where the gene is not detected at all is left empty). A black
     ring marks a coarse-family Fisher enrichment at BH FDR < 0.05
     (results/tables/geo_gse211630_candidate_celltype_enrichment.csv).
     Row labels are coloured by cell-type family.
  B  expression of each candidate on the OLP UMAP embedding. Non-expressing
     cells are drawn first in light grey, expressing cells on top in the same
     magma scale as panel A (0-4 log1p) so colour is comparable between panels;
     the enlarged marker is a visibility aid only and encodes no quantity.
  C  per-gene summary: fraction of cells expressing and mean log1p in the OLP
     atlas, in NEOLP and in EOLP, plus the maximum raw UMI observed

Inputs
  data/processed/geo_GSE211630/scanpy/gse211630_olp_lognorm.h5ad
  results/tables/geo_gse211630_celltype_annotation.csv
  results/tables/geo_gse211630_candidate_celltype_enrichment.csv
  results/tables/geo_gse211630_candidate_loc.csv
Output
  figures/supplementary/supp_fig08_candidate_expression.{png,pdf}  (300 dpi / vector)
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

exec(compile(open(r"\\?\D:\wisp-science\skills\figure-style\runtime.py",
                  encoding="utf-8").read(), "runtime.py", "exec"))

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
os.chdir(ROOT)
SCAN = os.path.join("data", "processed", "geo_GSE211630", "scanpy")
TAB = os.path.join("results", "tables")
FIGDIR = os.path.join("figures", "supplementary")
os.makedirs(FIGDIR, exist_ok=True)

apply_figure_style(sizes=(7, 6.3, 5.7))

CAND = ["GDF11", "SUOX", "EGFL8", "TAPBP", "HLA-DPA1"]
FAMILY_ORDER = ["T/NK", "B/Plasma", "Myeloid", "Mast", "Epithelial",
                "Fibroblast", "Pericyte/SMC", "Endothelial"]
FAM_COL = {"T/NK": "#C8102E", "B/Plasma": "#7B3FA0", "Myeloid": "#2E8B57",
           "Mast": "#B8860B", "Epithelial": "#E08214", "Fibroblast": "#3B6FB6",
           "Pericyte/SMC": "#5B3A29", "Endothelial": "#008B8B"}
# cluster -> short row label (cross-references the annotation table's Leiden ids;
# c12 and c18 are doublet clusters and carry no annotation)
SHORT = {0: "c0 Fibroblast SFRP2+", 1: "c1 Fibroblast CXCL14+", 2: "c2 Pericyte/SMC",
         3: "c3 T resting/recirculating", 4: "c4 T/NK stressed (HSP+)",
         5: "c5 T/NK low-complexity", 6: "c6 CD8+ T GZMK+ EM",
         7: "c7 CD4+ Treg", 8: "c8 T naive-like (IL7R+)",
         9: "c9 Plasmablast (GZMB+)", 10: "c10 Macrophage (LYZ+)",
         11: "c11 Epithelial basal/prolif.", 13: "c13 Mast",
         14: "c14 B cell", 15: "c15 Plasma cell",
         16: "c16 CD8+ T GZMB+/CXCL13+", 17: "c17 Endo lymphatic (CCL21+)",
         19: "c19 Endo blood (ENG+)", 20: "c20 Myeloid S100-hi"}
DOUBLET = (12, 18)
VMAX = 4.0
# panel A dot area in pt^2 at fraction = 1. Row pitch is ~11 pt, so a full
# marker (diameter ~7.7 pt) still overlaps nothing adjacent. S_MIN keeps a
# fraction detected in only a handful of cells visible (0.5 pt^2 ~ 0.011), so
# that an empty cell really does mean "no transcript detected".
S_MAX = 46.0
S_MIN = 0.5

# Original plate geometry (fractions of the 180 x 178 mm figure). main() reuses these so the
# published plate is reproduced exactly; the recomposition builder passes its own.
PLATE_W_MM, PLATE_H_MM = 180.0, 178.0
RECT_A = [0.315, 0.505, 0.395, 0.415]
RECT_A_CBAR = [0.755, 0.640, 0.020, 0.180]
RECT_A_LEGEND = [0.735, 0.455, 0.250, 0.168]
RECT_B_GRID = dict(left=0.045, right=0.895, bottom=0.235, top=0.445, wspace=0.06)
RECT_B_CBAR = [0.925, 0.255, 0.016, 0.150]
RECT_C = [0.045, 0.028, 0.930, 0.155]

# Panel titles (drawn as artist text in the original plate; kept here so a recomposed panel can
# suppress or restyle them).
TITLE_A = ("candidate expression across 19 OLP cell types "
           "(c12/c18 doublet clusters excluded)")

_data = None


def _ensure():
    """Load the atlas and derive panel inputs once per process.

    The h5ad read dominates runtime, so it is cached: the recomposition builder draws the same
    panels into several plates.
    """
    global _data
    if _data is not None:
        return _data

    import anndata as ad

    olp = ad.read_h5ad(os.path.join(SCAN, "gse211630_olp_lognorm.h5ad"))
    ann = pd.read_csv(os.path.join(TAB, "geo_gse211630_celltype_annotation.csv"))
    fam_of = {int(r["leiden_0.5"]): r["celltype_family"] for _, r in ann.iterrows()
              if isinstance(r["celltype_family"], str) and r["celltype_family"].strip()}
    enr = pd.read_csv(os.path.join(TAB, "geo_gse211630_candidate_celltype_enrichment.csv"))
    loc = pd.read_csv(os.path.join(TAB, "geo_gse211630_candidate_loc.csv")).set_index("gene")

    leid = olp.obs["leiden_0.5"].astype(int).values
    umap = olp.obsm["X_umap"] if "X_umap" in olp.obsm else \
        olp.obsm[[k for k in olp.obsm.keys() if "umap" in k.lower()][0]]

    missing = [g for g in CAND if g not in olp.var_names]
    print("candidates missing from the OLP object:", missing, flush=True)
    Xc = olp[:, CAND].X
    Xc = np.asarray(Xc.todense()) if hasattr(Xc, "todense") else np.asarray(Xc)

    cl_list = [c for c in sorted(fam_of)]
    assert len(cl_list) == 19, f"expected 19 annotated clusters, got {len(cl_list)}"
    rows = []
    for cl in cl_list:
        m = leid == cl
        for j, g in enumerate(CAND):
            v = Xc[m, j]
            ne = int((v > 0).sum())
            rows.append({"cluster": cl, "family": fam_of[cl], "gene": g,
                         "frac": ne / int(m.sum()),
                         "mean_expr": float(v[v > 0].mean()) if ne else 0.0})
    dp = pd.DataFrame(rows)
    sig = {(g, f) for g, f in zip(enr["gene"][enr["FDR" if "FDR" in enr else
                                             "fisher_FDR_BH"] < 0.05],
                                  enr["family"][enr["FDR" if "FDR" in enr else
                                                   "fisher_FDR_BH"] < 0.05])}
    print("significant (gene, family) pairs:", len(sig), flush=True)

    _data = dict(fam_of=fam_of, leid=leid, umap=umap, Xc=Xc, cl_list=cl_list,
                 dp=dp, sig=sig, loc=loc, missing=missing)
    return _data


def panel_A(fig, rect=None, cbar_rect=None, legend_rect=None, title=TITLE_A,
            with_letter=True):
    """Candidate x cell-type dot plot, plus its size legend and colour bar.

    Returns the axes it created so a caller can run a collision scan over them.
    """
    D = _ensure()
    dp, sig, fam_of = D["dp"], D["sig"], D["fam_of"]
    cl_list = D["cl_list"]
    rect = RECT_A if rect is None else list(rect)

    axA = fig.add_axes(rect)
    order = [c for f in FAMILY_ORDER for c in cl_list if fam_of[c] == f]
    ypos = {cl: i for i, cl in enumerate(reversed(order))}
    prev = None
    for cl in reversed(order):
        f = fam_of[cl]
        if prev is not None and f != prev:
            axA.axhline(ypos[cl] + 0.5, color="#CFCFCF", lw=0.6)
        prev = f
    for r in dp.itertuples():
        if float(r.frac) <= 0:
            continue  # not detected in this cell type: leave the cell empty
        ring = (r.gene, r.family) in sig
        axA.scatter(CAND.index(r.gene), ypos[r.cluster],
                    s=max(S_MIN, S_MAX * float(r.frac)),
                    c=[r.mean_expr], cmap="magma", vmin=0, vmax=VMAX,
                    linewidths=0.7 if ring else 0,
                    edgecolors="#111111" if ring else "none", alpha=0.94)
    axA.set_xticks(range(len(CAND)))
    axA.set_xticklabels(CAND, fontsize=6.5, style="italic")
    axA.set_yticks([ypos[c] for c in order])
    axA.set_yticklabels([SHORT[c] for c in order], fontsize=5.4)
    for t, c in zip(axA.get_yticklabels(), order):
        t.set_color(FAM_COL[fam_of[c]])
    axA.tick_params(axis="y", length=0, pad=1.5)
    axA.tick_params(axis="x", length=2, pad=1.5)
    axA.set_ylim(-0.6, 18.6)
    axA.set_xlim(-0.6, len(CAND) - 0.4)
    axA.set_ylabel("annotated cell type (Leiden cluster)", fontsize=5.7)
    set_frame(axA)
    if with_letter:
        panel_letter(axA, "A", dx=-0.62, dy=1.045, fontsize=8.5)
    if title:
        axA.text(0.0, 1.012, title,
                 transform=axA.transAxes, fontsize=5.4, color="#444444", va="bottom")

    cb = cax = None
    if cbar_rect is not None:
        cax = fig.add_axes(list(cbar_rect))
        cb = fig.colorbar(plt.cm.ScalarMappable(
            norm=matplotlib.colors.Normalize(0, VMAX), cmap="magma"), cax=cax)
        cb.set_label("mean log1p expr.\n(expressing cells)", fontsize=5.4, labelpad=2.0)
        cb.ax.tick_params(labelsize=5.2, length=2, pad=1.5)

    axL = None
    if legend_rect is not None:
        # matplotlib's legend spaces handles by handleheight x fontsize, which cannot
        # hold markers several times that tall: the size key is drawn on its own axes.
        axL = fig.add_axes(list(legend_rect))
        axL.axis("off")
        axL.set_xlim(0, 1)
        axL.set_ylim(0, 1)
        axL.set_xticks([])
        axL.set_yticks([])
        axL.text(0.0, 0.995, "frac. expressing", fontsize=5.4, va="top")
        for i, (lab, area) in enumerate([("1", S_MAX), ("0.75", S_MAX * 0.75),
                                         ("0.5", S_MAX * 0.5), ("0.25", S_MAX * 0.25),
                                         ("\u22640.01", S_MIN)]):
            yv = 0.80 - i * 0.155
            axL.scatter([0.14], [yv], s=area, c="#8C8C8C", linewidths=0)
            axL.text(0.34, yv, lab, fontsize=5.4, va="center", ha="left")
        axL.scatter([0.14], [0.03], s=S_MAX * 0.55, facecolors="none",
                    edgecolors="#111111", linewidths=0.7)
        axL.text(0.34, 0.03, "Fisher FDR<0.05", fontsize=5.4, va="center", ha="left")

    return dict(axA=axA, cax=cax, cb=cb, axL=axL)


def panel_B(fig, rect=None, cbar_rect=None, with_letter=True, ncol=5):
    """Per-candidate expression on the OLP UMAP, one facet per gene."""
    D = _ensure()
    Xc, umap, loc = D["Xc"], D["umap"], D["loc"]
    if rect is None:
        gsB = fig.add_gridspec(1, 5, **RECT_B_GRID)
    else:
        gsB = fig.add_gridspec(1, ncol, left=rect[0], right=rect[0] + rect[2],
                               bottom=rect[1], top=rect[1] + rect[3], wspace=0.06)
    normB = matplotlib.colors.Normalize(0, VMAX)
    cmapB = plt.get_cmap("magma")
    b_ax = []
    for j, g in enumerate(CAND):
        ax = fig.add_subplot(gsB[0, j])
        b_ax.append(ax)
        v = Xc[:, j]
        on = v > 0
        # background first, expressing cells on top: one scatter would let the later
        # non-expressing cells in the array order overdraw a rare gene's signal
        ax.scatter(umap[~on, 0], umap[~on, 1], s=0.06, c="#F0F0F0", linewidths=0,
                   rasterized=True)
        ax.scatter(umap[on, 0], umap[on, 1], s=0.5, c=v[on], cmap=cmapB, norm=normB,
                   linewidths=0, rasterized=True)
        fr = float(loc.loc[g, "frac_cells"])
        ax.set_title(f"{g}\n{100 * fr:.1f}% of cells", fontsize=5.7, style="italic")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_box_aspect(0.95)
        set_frame(ax)
        if j == 0:
            ax.set_ylabel("UMAP 2", fontsize=5.4)
        if j == 2:
            ax.set_xlabel("UMAP 1", fontsize=5.4)
    if with_letter:
        panel_letter(b_ax[0], "B", dx=-0.14, dy=1.22, fontsize=8.5)

    cbB = caxB = None
    if cbar_rect is not None:
        caxB = fig.add_axes(list(cbar_rect))
        cbB = fig.colorbar(plt.cm.ScalarMappable(norm=normB, cmap=cmapB), cax=caxB)
        cbB.set_label("log1p expr.\nper cell", fontsize=5.4, labelpad=2.0)
        cbB.ax.tick_params(labelsize=5.2, length=2, pad=1.5)

    return dict(b_ax=b_ax, caxB=caxB, cbB=cbB)


def panel_C(fig, rect=None, with_letter=True):
    """Per-gene detection/expression summary table."""
    D = _ensure()
    loc = D["loc"]
    axC = fig.add_axes(RECT_C if rect is None else list(rect))
    axC.axis("off")
    # an axis-free axes autoscales its data limits to the axes-fraction artists drawn
    # into it, which would push the transData table text off the canvas: pin the unit
    # square and place every text in axes fraction.
    axC.set_xlim(0, 1)
    axC.set_ylim(0, 1)
    axC.set_xticks([])
    axC.set_yticks([])
    cols = ["gene", "% cells\nOLP", "% cells\nNEOLP", "% cells\nEOLP",
            "mean log1p\nNEOLP", "mean log1p\nEOLP", "max raw\nUMI"]
    xs = [0.02, 0.20, 0.33, 0.46, 0.59, 0.74, 0.89]
    for x, c in zip(xs, cols):
        axC.text(x, 0.93, c, fontsize=5.4, va="center", ha="left", color="#222222",
                 transform=axC.transAxes)
    axC.plot([0.0, 1.0], [0.855, 0.855], color="#BBBBBB", lw=0.6,
             transform=axC.transAxes, clip_on=False)
    for i, g in enumerate(CAND):
        y = 0.74 - i * 0.145
        vals = [g,
                f"{100 * loc.loc[g, 'frac_cells']:.1f}",
                f"{100 * loc.loc[g, 'NEOLP_frac']:.1f}",
                f"{100 * loc.loc[g, 'EOLP_frac']:.1f}",
                f"{loc.loc[g, 'NEOLP_mean_log1p']:.3f}",
                f"{loc.loc[g, 'EOLP_mean_log1p']:.3f}",
                f"{int(loc.loc[g, 'raw_max_umi_olp'])}"]
        for k, (x, v) in enumerate(zip(xs, vals)):
            axC.text(x, y, v, fontsize=5.4, va="center", ha="left",
                     style="italic" if k == 0 else "normal",
                     transform=axC.transAxes)
        if i < len(CAND) - 1:
            axC.plot([0.0, 1.0], [y - 0.072, y - 0.072], color="#EEEEEE", lw=0.5,
                     transform=axC.transAxes, clip_on=False)
    if with_letter:
        panel_letter(axC, "C", dx=0.0, dy=1.03, fontsize=8.5)
    return dict(axC=axC)


def collision_scan(fig, A, B, C):
    """Assert that no text overlaps, leaves the canvas, or crosses into another panel."""
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()

    axA, cax, cb, axL = A["axA"], A["cax"], A["cb"], A["axL"]
    b_ax, caxB, cbB = B["b_ax"], B["caxB"], B["cbB"]
    axC = C["axC"]

    def texts_of(ax):
        out = [t for t in ax.texts if t.get_text().strip()]
        for t in (ax.title, ax.xaxis.label, ax.yaxis.label):
            if t.get_text().strip():
                out.append(t)
        # an axis that is switched off still owns tick label objects: scanning them
        # would report collisions between text that is never drawn.
        if ax.xaxis.get_visible() and ax.get_xticks().size:
            out += [t for t in ax.get_xticklabels() if t.get_text().strip()]
        if ax.yaxis.get_visible() and ax.get_yticks().size:
            out += [t for t in ax.get_yticklabels() if t.get_text().strip()]
        return out

    all_ax = [axA] + b_ax + [axC] + [a for a in (cax, axL, caxB) if a is not None]
    arts_raw = []
    for a_ in all_ax:
        arts_raw += texts_of(a_)
    for cab in (cb, cbB):
        if cab is None:
            continue
        if cab.ax.yaxis.label.get_text().strip():
            arts_raw.append(cab.ax.yaxis.label)
        arts_raw += [t for t in cab.ax.get_yticklabels() if t.get_text().strip()]
    for leg in fig.legends:
        arts_raw += [t for t in leg.get_texts() if t.get_text().strip()]
        if leg.get_title().get_text().strip():
            arts_raw.append(leg.get_title())
    # a colourbar axes may be reached twice, so collapse repeats by artist identity
    # before the pairwise comparison.
    _seen = set()
    arts = []
    for t in arts_raw:
        if id(t) in _seen:
            continue
        _seen.add(id(t))
        arts.append(t)

    TOL = 2.5

    def ov(ba, bb):
        return (min(ba.x1, bb.x1) - max(ba.x0, bb.x0) > TOL and
                min(ba.y1, bb.y1) - max(ba.y0, bb.y0) > TOL)

    hits = []
    boxes = [(t, t.get_window_extent(renderer=rend)) for t in arts]
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            if ov(boxes[i][1], boxes[j][1]):
                hits.append((boxes[i][0].get_text(), boxes[j][0].get_text()))
    fb = fig.get_window_extent(renderer=rend)
    outc = [t.get_text() for t, b in boxes
            if b.x0 < -TOL or b.y0 < -TOL or b.x1 > fb.x1 + TOL or b.y1 > fb.y1 + TOL]

    # a panel's own text must not run into a different panel: text/text scanning
    # alone cannot see a label crossing an axes it does not belong to.
    others_C = [axC] + ([axL] if axL is not None else [])
    others_B = b_ax + [axC] + ([axL] if axL is not None else [])
    intr = []
    for owner, others in ((axA, b_ax + others_C), (axC, [axA] + b_ax),
                          (axL, [axA] + b_ax + [axC]),
                          (b_ax[0], [axC]), (b_ax[-1], [axC] + ([caxB] if caxB is not None else []))):
        if owner is None:
            continue
        for t in texts_of(owner):
            b = t.get_window_extent(renderer=rend)
            for oa in others:
                if ov(b, oa.get_window_extent(renderer=rend)):
                    intr.append((t.get_text(), str(oa.get_position().bounds)))
    print("text objects scanned:", len(arts), flush=True)
    print("overlapping text pairs:", hits, flush=True)
    print("text outside canvas:", outc, flush=True)
    print("labels crossing into another panel:", intr, flush=True)
    assert not hits, "collision scan failed (text/text)"
    assert not outc, "collision scan failed (off-canvas text)"
    assert not intr, "collision scan failed (label vs neighbouring panel)"


def main():
    MM = 1 / 25.4
    fig = plt.figure(figsize=(PLATE_W_MM * MM, PLATE_H_MM * MM))

    A = panel_A(fig, RECT_A, RECT_A_CBAR, RECT_A_LEGEND)
    B = panel_B(fig, None, RECT_B_CBAR)
    C = panel_C(fig, RECT_C)

    collision_scan(fig, A, B, C)

    fig.savefig(os.path.join(FIGDIR, "supp_fig08_candidate_expression.png"), dpi=300,
                bbox_inches="tight")
    fig.savefig(os.path.join(FIGDIR, "supp_fig08_candidate_expression.pdf"),
                bbox_inches="tight")
    D = _ensure()
    order = [c for f in FAMILY_ORDER for c in D["cl_list"] if D["fam_of"][c] == f]
    print("rows in A:", len(order), "| significant pairs:", len(D["sig"]), flush=True)
    print("frac range in A:", round(D["dp"]["frac"].min(), 4),
          round(D["dp"]["frac"].max(), 4), flush=True)
    print("saved supp_fig08_candidate_expression.{png,pdf}", flush=True)


if __name__ == "__main__":
    main()
