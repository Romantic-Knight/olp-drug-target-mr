# -*- coding: utf-8 -*-
"""Supplementary Figure 6: GSE211630 OLP atlas — raw Leiden clustering and
per-sample coverage on the analysed 46,478-cell OLP embedding.

Why these two views live in the Supplement (not the main text)
---------------------------------------------------------------
* Main-text Figure 6A shows the *annotated* atlas (19 subtypes, family
  palette) and Figure 6B the NEOLP-vs-EOLP condition overlay.
* This figure keeps the two remaining overview views that are useful context
  but not load-bearing for the erosion/composition/communication analyses:
  (A) the un-annotated Leiden res-0.5 clusters c0-c20 (cluster colours are
  arbitrary and independent of the family palette; doublet-suspect c12/c18
  grey, n=32, excluded from all statistics), and (B) per-sample cell coverage
  on the same embedding (all 46,478 sequenced OLP cells; NEOLP-1/2/3 then
  EOLP-1/2), illustrating that neither condition nor any single sample drives
  a cluster.

Outputs
-------
figures/supplementary/supp_fig06_olp_umap.png (300 dpi) + .pdf
analysis/scrna/_qa/supp_fig06_olp_umap.png (150 dpi QA copy)

Panel API (2026-09-23)
----------------------
Refactored for figure recombination (Plan 3): each panel is now a function that draws onto an
Axes supplied by the caller, so a composing script can place the same panel in a new plate. The
drawing code is unchanged; `main()` reproduces the original two-panel figure byte-for-byte, so the
canonical output stays reproducible.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FIG = os.path.join(ROOT, "figures", "supplementary")
QA = os.path.join(ROOT, "analysis", "scrna", "_qa")
H5 = os.path.join(ROOT, "data", "processed", "geo_GSE211630", "scanpy",
                  "gse211630_olp_lognorm.h5ad")
OB = os.path.join(ROOT, "data", "processed", "geo_GSE211630", "scanpy",
                  "gse211630_olp_obs_annotation.csv")

SAMP_ORDER = ["NEOLP-1", "NEOLP-2", "NEOLP-3", "EOLP-1", "EOLP-2"]
SAMP_COL = {"NEOLP-1": "#8FB8DE", "NEOLP-2": "#4477AA", "NEOLP-3": "#2A5684",
            "EOLP-1": "#F7B285", "EOLP-2": "#EE7733"}

_DATA = None


def load():
    """Load the OLP embedding and annotation once per process."""
    global _DATA
    if _DATA is not None:
        return _DATA
    olp = sc.read_h5ad(H5)                       # 46,478 x 17,669, X_umap ready
    ob = pd.read_csv(OB, index_col=0)
    assert (olp.obs.index == ob.index).all()
    olp.obs["celltype_sub"] = ob["celltype_sub"].values
    umap = np.asarray(olp.obsm["X_umap"])
    leid = olp.obs["leiden_0.5"].astype(int).values
    sample = olp.obs["sample"].astype(str).values
    assert umap.shape[0] == 46478 and olp.obs["sample"].nunique() == 5

    cl_sub = {}
    for cl in sorted(set(leid) - {12, 18}):
        m = leid == cl
        cl_sub[cl] = str(olp.obs["celltype_sub"][m].dropna().iloc[0])
    lab = {cl: f"c{cl}  {cl_sub[cl]}" for cl in sorted(cl_sub)}
    lab.update({12: "c12  doublet-suspect (excluded)",
                18: "c18  doublet-suspect (excluded)"})

    # arbitrary per-cluster colours (independent of the family palette of Fig 8A)
    import colorsys
    cmap = {}
    k = 0
    for cl in range(21):
        if cl in (12, 18):
            cmap[cl] = (0.72, 0.72, 0.72)
        else:
            h = (0.02 + k * 0.618034) % 1.0     # golden-angle hue sweep, no greys
            cmap[cl] = colorsys.hls_to_rgb(h, 0.52, 0.66)
            k += 1
    _DATA = {
        "umap": umap, "leid": leid, "sample": sample, "lab": lab, "cmap": cmap,
        "gry": (0.72, 0.72, 0.72),
        "samp_n": {s: int((sample == s).sum()) for s in SAMP_ORDER},
    }
    return _DATA


def strip_axes(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def panel_a(ax, tag=True, letter="A", ncol=1, labelspacing=0.26, title=None,
            ms=6.0, ms_dbl=5.2, legend_fontsize=5.8):
    """(A) raw Leiden clusters c0-c20, doublet-suspect clusters in grey."""
    D = load()
    umap, leid, cmap, gry = D["umap"], D["leid"], D["cmap"], D["gry"]
    strip_axes(ax)
    dbl = np.isin(leid, (12, 18))
    for cl in range(21):
        m = (leid == cl) & ~dbl
        if m.sum() == 0:
            continue
        ax.scatter(umap[m, 0], umap[m, 1], s=1.5, c=[cmap[cl]] * int(m.sum()),
                   linewidths=0, alpha=0.85, rasterized=True)
    ax.scatter(umap[dbl, 0], umap[dbl, 1], s=3.8, c=[gry] * int(dbl.sum()),
               linewidths=0.4, edgecolors="#444444", alpha=0.98, rasterized=True)
    ax.set_title(title or "Raw Leiden clustering (resolution 0.5, clusters c0–c20)",
                 loc="left", fontsize=8.2, weight="bold", pad=3)
    handles, labels = [], []
    for cl in sorted(D["lab"]):
        if cl in (12, 18):
            continue
        handles.append(plt.Line2D([], [], marker="o", ls="", markersize=ms,
                                  color=cmap[cl], alpha=0.95))
        labels.append(D["lab"][cl])
    handles.append(plt.Line2D([], [], marker="o", ls="", markersize=ms_dbl, color=gry,
                              markeredgecolor="#444444", markeredgewidth=0.6))
    labels.append("c12/c18 doublet-suspect (n = 32, excluded)")
    ax.legend(handles, labels, loc="upper left", bbox_to_anchor=(1.02, 1.0),
              frameon=True, framealpha=0.96, edgecolor="#DDDDDD",
              fontsize=legend_fontsize,
              labelspacing=labelspacing, handletextpad=0.45, borderpad=0.5, ncol=ncol,
              columnspacing=1.0)
    if tag:
        ax.text(-0.030, 1.05, letter, transform=ax.transAxes, fontsize=11,
                fontweight="bold", ha="left", va="bottom")
    return ax


def panel_b(ax, tag=True, letter="B", title=None, legend_loc="lower left",
            legend_fontsize=6.0, ms=5.6, labelspacing=None, borderpad=None,
            legend_bbox=None, legend_ncol=1):
    """(B) per-sample coverage over all 46,478 sequenced OLP cells."""
    D = load()
    umap, sample = D["umap"], D["sample"]
    strip_axes(ax)
    for s in SAMP_ORDER:
        m = sample == s
        ax.scatter(umap[m, 0], umap[m, 1], s=1.5, c=[SAMP_COL[s]] * int(m.sum()),
                   linewidths=0, alpha=0.85, rasterized=True)
    ax.set_title(title or "Per-sample coverage (all sequenced OLP cells)",
                 loc="left", fontsize=8.2, weight="bold", pad=3)
    hb = [plt.Line2D([], [], marker="o", ls="", markersize=ms, color=SAMP_COL[s],
                     label=f"{s} (n = {D['samp_n'][s]:,})") for s in SAMP_ORDER]
    # legend_loc is a parameter because the map's busy quadrant differs between plates: the original
    # standalone panel put the key at lower-left, but in a composite that position covers the
    # embedding's lower-left structure.
    # labelspacing / borderpad default to the source values; a denser recomposition can pass
    # smaller ones to keep the key from growing across the embedding.
    # legend_bbox / legend_ncol let a recomposition park the key outside the axes (e.g. in the band
    # under the panel) when no corner of the embedding is empty enough to hold it.
    ax.legend(handles=hb, loc=legend_loc, frameon=True, framealpha=0.9,
              edgecolor="#DDDDDD", fontsize=legend_fontsize, handletextpad=0.4,
              labelspacing=0.5 if labelspacing is None else labelspacing,
              borderpad=0.4 if borderpad is None else borderpad,
              bbox_to_anchor=legend_bbox, ncol=legend_ncol,
              columnspacing=0.9 if legend_ncol > 1 else 2.0)
    if tag:
        ax.text(-0.030, 1.05, letter, transform=ax.transAxes, fontsize=11,
                fontweight="bold", ha="left", va="bottom")
    return ax


def build(width_in=13.0, height_in=6.4):
    """Original two-panel plate (square-ish boxes, panels left and right)."""
    fig = plt.figure(figsize=(width_in, height_in))
    S = 3.95                                    # panel side in inches
    yC = 0.475
    panel_a(fig.add_axes([0.012, yC - (S / height_in) / 2, S / width_in, S / height_in]))
    panel_b(fig.add_axes([0.635, yC - (S / height_in) / 2, S / width_in, S / height_in]))
    return fig


def main():
    os.makedirs(QA, exist_ok=True)
    fig = build()
    fig.savefig(os.path.join(FIG, "supp_fig06_olp_umap.png"), dpi=300)
    fig.savefig(os.path.join(FIG, "supp_fig06_olp_umap.pdf"))
    fig.savefig(os.path.join(QA, "supp_fig06_olp_umap.png"), dpi=150)
    plt.close(fig)
    print("supp fig06 written")


if __name__ == "__main__":
    main()
