# -*- coding: utf-8 -*-
"""Supplementary Figure 5 - GSE211630 dataset provenance and quality control.

New plate built when the supplementary set was reordered to follow the
manuscript narrative (2026-09-13): dataset provenance and QC now precede the
atlas clustering view, which became Supplementary Figure 6.

Panels
  A  per-sample yield in the raw matrices: cells, detected gene rows and mean
     UMI per cell (6 samples: 1 healthy control, 3 NEOLP, 2 EOLP)
  B  per-cell QC per sample (log10 genes, log10 counts and the sample-specific
     filter metric actually applied to that sample - Scrublet doublet score
     where doublet detection was run, mitochondrial fraction otherwise),
     kept versus excluded cells
  C  per-sample coverage on the all-sample UMAP embedding (60,308 cells)

Inputs
  results/tables/geo_gse211630_qc_summary.csv          (panel A)
  data/processed/geo_GSE211630/qc_cell_metrics.csv.gz  (panel B, 62,357 cells)
  data/processed/geo_GSE211630/preprocess_summary.json (panels A and B)
  data/processed/geo_GSE211630/scanpy/gse211630_all_analyzed.h5ad (panel C)

Output
  figures/supplementary/supp_fig05_dataset_qc.{png,pdf}   (300 dpi / vector)

Panel API (2026-09-23)
----------------------
Refactored for figure recombination (Plan 3). Each panel is now a function that fills a
caller-supplied gridspec cell, creating its own subgridspec inside it:

    panel_a(fig, cell)   -> [3 axes];  legend_a(fig, cell) draws its colour key
    panel_b(fig, cell)   -> [18 axes]; legend_b(fig, cell) draws its kept/excluded key
    panel_c(fig, cell)   -> [6 axes]

Because every panel builds its axes from the cell it is handed, a composing script can give the
same panel a different area without touching the drawing code. `main()` reproduces the original
three-panel plate, verified byte-identical to the pre-refactor output.
"""
import json
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
OUT = os.path.join("data", "processed", "geo_GSE211630")
SCAN = os.path.join(OUT, "scanpy")
TAB = os.path.join("results", "tables")
FIGDIR = os.path.join("figures", "supplementary")
os.makedirs(FIGDIR, exist_ok=True)

apply_figure_style(sizes=(7, 6.3, 5.7))

# ---------------------------------------------------------------- house colours
C_NORMAL, C_NEOLP, C_EOLP = "#9E9E9E", "#3B6FB6", "#E08214"
GRP_COL = {"normal": C_NORMAL, "NEOLP": C_NEOLP, "EOLP": C_EOLP}
GRP_NAME = {"normal": "normal (n = 1)", "NEOLP": "NEOLP (n = 3)",
            "EOLP": "EOLP (n = 2)"}
C_KEPT, C_EXCL = "#4D4D4D", "#C8102E"
ORDER = ["normal-1", "NEOLP-1", "NEOLP-2", "NEOLP-3", "EOLP-1", "EOLP-2"]

# ---------------------------------------------------------------- load inputs
qc = pd.read_csv(os.path.join(TAB, "geo_gse211630_qc_summary.csv"))
qc = qc.set_index("label").loc[ORDER].reset_index()

cell = pd.read_csv(os.path.join(OUT, "qc_cell_metrics.csv.gz"))
with open(os.path.join(OUT, "preprocess_summary.json"), encoding="utf-8") as f:
    pre = json.load(f)
MT_LABELS = {v["label"] for g, v in pre["samples"].items()
             if g in set(pre["params"]["mt_samples"])}
DOUBLET_THR = {v["label"]: (None if v["doublet_threshold"] is None
                            else float(v["doublet_threshold"]))
               for v in pre["samples"].values()}
TOT_RAW = pre["totals"]["total_cells_raw"]
TOT_KEPT = pre["totals"]["total_cells_kept"]

import anndata as ad
allobj = ad.read_h5ad(os.path.join(SCAN, "gse211630_all_analyzed.h5ad"))

rng = np.random.default_rng(1337)


def sub(n, k):
    return rng.choice(n, size=min(k, n), replace=False)


def violin(ax, v):
    ax.violinplot(v, positions=[0], widths=0.86, showextrema=False)
    for pc in ax.collections:
        pc.set_facecolor("#D9D9D9")
        pc.set_edgecolor("#8C8C8C")
        pc.set_linewidth(0.35)
        pc.set_alpha(0.95)


def dots(ax, v, kept):
    """kept cells as a grey strip sample, excluded cells as red points."""
    ix = sub(len(v), 600)
    ax.scatter(rng.uniform(-0.11, 0.11, len(ix)), v[ix], s=0.5, c=C_KEPT,
               alpha=0.35, linewidths=0, rasterized=True)
    j = np.flatnonzero(~kept)
    ax.scatter(rng.uniform(-0.20, 0.20, len(j)), v[j], s=1.1, c=C_EXCL,
               alpha=0.55, linewidths=0, rasterized=True)


# ================================================================ A: per sample
def panel_a(fig, cell_gs, with_letter=True):
    """(A) per-sample yield: cells, gene rows and mean UMI. Returns the three axes."""
    gsA = cell_gs.subgridspec(1, 3, wspace=0.42)
    a_ax = []
    for k, (col, title, ylab) in enumerate([
            ("n_cells", "Cells per sample", "cells"),
            ("n_genes", "Gene rows per sample", "gene rows"),
            ("mean_umi_per_cell", "Mean UMI per cell", "UMI")]):
        ax = fig.add_subplot(gsA[0, k])
        a_ax.append(ax)
        xy = np.arange(len(qc))
        ax.bar(xy, qc[col].values, width=0.72,
               color=[GRP_COL[g] for g in qc["group"]], edgecolor="none")
        ax.set_xticks(xy)
        ax.set_xticklabels(qc["label"], rotation=45, ha="right", fontsize=5.0)
        ax.set_title(title)
        ax.set_ylabel(ylab)
        top = float(qc[col].max())
        ax.set_ylim(0, top * 1.26)
        for x, v in zip(xy, qc[col].values):
            ax.text(x, v + top * 0.04, f"{v:,.0f}", ha="center", va="bottom",
                    fontsize=5.0, color="#222222")
        set_frame(ax)
    if with_letter:
        panel_letter(a_ax[0], "A", dx=-0.26, dy=1.18, fontsize=8.5)
    return a_ax


def legend_a(fig, cell_gs):
    """Colour key for panel A (sits in its own thin cell of the plate)."""
    handles = [plt.Rectangle((0, 0), 1, 1, fc=GRP_COL[g], ec="none")
               for g in ("normal", "NEOLP", "EOLP")]
    cellA = cell_gs.get_position(fig)
    fig.legend(handles, [GRP_NAME[g] for g in ("normal", "NEOLP", "EOLP")],
               loc="center", bbox_to_anchor=(0.5, cellA.y0 + cellA.height / 2.0),
               ncol=3, handlelength=1.0, handletextpad=0.5, columnspacing=1.6)


# ================================================================ B: per cell
def panel_b(fig, cell_gs, with_letter=True):
    """(B) per-cell QC per sample. Returns the 18 axes (3 rows x 6 samples)."""
    gsB = cell_gs.subgridspec(3, 6, wspace=0.30, hspace=0.58)
    b_ax = []
    for r in range(3):
        for c in range(6):
            b_ax.append(fig.add_subplot(gsB[r, c]))
    by_lab = cell.groupby("label")
    for c, s in enumerate(ORDER):
        d = by_lab.get_group(s)
        for r, col in enumerate(["n_genes", "n_counts"]):
            ax = b_ax[r * 6 + c]
            violin(ax, np.log10(d[col].values.astype(float) + 1.0))
            dots(ax, np.log10(d[col].values.astype(float) + 1.0),
                 d["kept_final"].values)
            ax.set_xlim(-0.55, 0.55)
            ax.set_xticks([])
            set_frame(ax)
            if c == 0:
                ax.set_ylabel(["log10 genes", "log10 counts"][r])
            else:
                ax.set_yticklabels([])
        # row 3: the filter metric actually applied to this sample
        ax = b_ax[2 * 6 + c]
        if s in MT_LABELS:
            v = d["mt_frac"].values.astype(float) * 100.0
            violin(ax, v)
            dots(ax, v, d["kept_final"].values)
            ax.axhline(20.0, color=C_EXCL, lw=0.7, ls="--")
            ax.axhline(10.0, color="#8C8C8C", lw=0.7, ls=":")
            ax.set_ylim(0, 56)
            ax.text(0.04, 0.99, f"excl >20%: {int((~d['kept_final'].values).sum())}"
                                f"\nref ≥10%: {int((v >= 10.0).sum())}",
                    fontsize=5.0, transform=ax.transAxes, va="top", ha="left",
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.82,
                              pad=1.0))
            ax.set_title("mito %")
        else:
            v = d["doublet_score"].values.astype(float)
            violin(ax, v)
            dots(ax, v, d["kept_final"].values)
            thr = DOUBLET_THR[s]
            ax.axhline(thr, color=C_EXCL, lw=0.7, ls="--")
            ax.text(0.04, 0.99, f"pred: {int(d['predicted_doublet'].sum())}"
                                f"\nthr: {thr:.2f}"
                                f"\nexcl: {int((~d['kept_final'].values).sum())}",
                    fontsize=5.0, transform=ax.transAxes, va="top", ha="left",
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.82,
                              pad=1.0))
            ax.set_ylim(-0.06, 1.28)
            ax.set_title("doublet score")
        ax.set_xlim(-0.55, 0.55)
        ax.set_xticks([])
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
        set_frame(ax)
        b_ax[c].set_title(f"{s}\nkept {int(d['kept_final'].sum()):,} / raw {len(d):,}",
                          fontsize=5.5)
    if with_letter:
        panel_letter(b_ax[0], "B", dx=-0.52, dy=1.62, fontsize=8.5)
    return b_ax


def legend_b(fig, cell_gs):
    """kept / excluded key for panel B (sits in its own thin cell of the plate)."""
    handles = [plt.Line2D([], [], marker="o", ls="", ms=2.6, mfc=C_KEPT, mec="none"),
               plt.Line2D([], [], marker="o", ls="", ms=2.6, mfc=C_EXCL, mec="none")]
    cellB = cell_gs.get_position(fig)
    fig.legend(handles, ["kept", "excluded"], loc="center",
               bbox_to_anchor=(0.5, cellB.y0 + cellB.height / 2.0), ncol=2,
               handletextpad=0.4, columnspacing=1.8)


# ================================================================ C: coverage
def panel_c(fig, cell_gs, with_letter=True):
    """(C) per-sample coverage on the all-sample UMAP. Returns the six axes."""
    um = allobj.obsm["X_umap"]
    slab = allobj.obs["label"].astype(str).values
    gsC = cell_gs.subgridspec(2, 3, wspace=0.06, hspace=0.14)
    c_ax = []
    for k, s in enumerate(ORDER):
        ax = fig.add_subplot(gsC[k // 3, k % 3])
        c_ax.append(ax)
        g = qc.loc[qc["label"] == s, "group"].iloc[0]
        m = slab == s
        ax.scatter(um[~m, 0], um[~m, 1], s=0.06, c="#F0F0F0", linewidths=0,
                   rasterized=True)
        ax.scatter(um[m, 0], um[m, 1], s=0.06, c=GRP_COL[g], linewidths=0,
                   rasterized=True)
        ax.set_title(f"{s} (n = {int(m.sum()):,})", fontsize=5.7)
        if k >= 3:
            ax.set_xlabel("UMAP 1")
        if k % 3 == 0:
            ax.set_ylabel("UMAP 2")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_box_aspect(0.86)
        set_frame(ax)
    if with_letter:
        panel_letter(c_ax[0], "C", dx=-0.10, dy=1.10, fontsize=8.5)
    return c_ax


def collision_scan(fig):
    """Text/text and off-canvas QA over every axes of the plate (unchanged checks)."""
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    arts = []
    for ax in fig.axes:
        for t in ax.texts:
            if t.get_text().strip():
                arts.append(t)
        for t in (ax.title, ax.xaxis.label, ax.yaxis.label):
            if t.get_text().strip():
                arts.append(t)
        if ax.get_xticks().size:
            arts += [t for t in ax.get_xticklabels() if t.get_text().strip()]
        if len(ax.get_yticks()):
            arts += [t for t in ax.get_yticklabels() if t.get_text().strip()]
    for leg in fig.legends:
        arts += [t for t in leg.get_texts() if t.get_text().strip()]
    TOL = 2.5
    hits = []
    boxes = [(t, t.get_window_extent(renderer=rend)) for t in arts]
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i][1], boxes[j][1]
            if (min(a.x1, b.x1) - max(a.x0, b.x0) > TOL and
                    min(a.y1, b.y1) - max(a.y0, b.y0) > TOL):
                hits.append((boxes[i][0].get_text(), boxes[j][0].get_text()))
    fb = fig.get_window_extent(renderer=rend)
    out = [t.get_text() for t, b in boxes
           if b.x0 < -TOL or b.y0 < -TOL or b.x1 > fb.x1 + TOL or b.y1 > fb.y1 + TOL]
    print("text objects scanned:", len(arts), flush=True)
    print("overlapping text pairs:", hits, flush=True)
    print("text outside canvas:", out, flush=True)
    assert not hits, "collision scan failed (text/text)"
    assert not out, "collision scan failed (off-canvas text)"


def build():
    """Original 180 x 252 mm plate: A + its key, B + its key, C (5-row canvas)."""
    MM = 1 / 25.4
    fig = plt.figure(figsize=(180 * MM, 252 * MM))
    gs = fig.add_gridspec(5, 1, height_ratios=[1.00, 0.15, 2.05, 0.15, 2.45],
                          hspace=0.42, left=0.075, right=0.985, top=0.925,
                          bottom=0.030)
    panel_a(fig, gs[0])
    legend_a(fig, gs[1])
    panel_b(fig, gs[2])
    legend_b(fig, gs[3])
    panel_c(fig, gs[4])
    return fig


def main():
    fig = build()
    collision_scan(fig)
    fig.savefig(os.path.join(FIGDIR, "supp_fig05_dataset_qc.png"), dpi=300,
                bbox_inches="tight")
    fig.savefig(os.path.join(FIGDIR, "supp_fig05_dataset_qc.pdf"),
                bbox_inches="tight")
    print("samples:", len(qc), "| total_raw", TOT_RAW, "total_kept", TOT_KEPT, flush=True)
    print("excluded per sample:",
          cell.groupby("label")["kept_final"].apply(lambda s: int((~s).sum())).to_dict(),
          flush=True)


if __name__ == "__main__":
    main()
