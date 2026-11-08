# -*- coding: utf-8 -*-
"""Cell-level QC diagnosis figure for GSE211630 (figures/geo_gse211630_cellqc.png).

Rows: log10 n_genes | log10 n_counts | applied per-sample filter metric
      (mito% for the 3 MT-row samples; scrublet doublet score for the 3 MT-free).
Columns: 6 samples (normal-1, NEOLP-1/2/3, EOLP-1/2). Kept cells grey violins,
excluded cells red points; filter cut-offs drawn per panel. Honest n.a. markers
for the filter-metric row where a metric is undefined.
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
os.chdir(ROOT)
OUT = os.path.join("data", "processed", "geo_GSE211630")
FIG = os.path.join("figures", "geo_gse211630_cellqc.png")

# figure-style sidecar (pure geometry + role-mapped font ladder)
sidecar = r"\\?\D:\wisp-science\skills\figure-style\kernel.py"
exec(compile(open(sidecar, encoding="utf-8").read(), "kernel.py", "exec"))
apply_figure_style(sizes=(8, 7, 6))

m = pd.read_csv(os.path.join(OUT, "qc_cell_metrics.csv.gz"))
with open(os.path.join(OUT, "preprocess_summary.json"), encoding="utf-8") as f:
    summ = json.load(f)
samples = list(summ["samples"].keys())
labels = [summ["samples"][g]["label"] for g in samples]
KEPT_C, EXCL_C = "#B9B9B9", "#B2182B"

# --- per-sample doublet-score diagnostics (for the delivery message) ---
dbl = m[m["sample"].isin(["GSM6481637", "GSM6481638", "GSM6481642"])]
print("doublet score quantiles (all cells):", flush=True)
print(dbl.groupby("sample")["doublet_score"].quantile([0.5, 0.9, 0.99, 1.0])
      .round(3).to_string(), flush=True)

rng = np.random.default_rng(7)
fig, axes = plt.subplots(3, 6, figsize=(12.6, 8.0), sharex=True)
fig.subplots_adjust(left=0.055, right=0.995, top=0.96, bottom=0.11,
                    wspace=0.16, hspace=0.34)

rows = [("n_genes", "log10 n genes", None),
        ("n_counts", "log10 n counts", None)]

for ci, gsm in enumerate(samples):
    sub = m[m["sample"] == gsm]
    kept = sub[sub.kept_final]
    excl = sub[~sub.kept_final]
    meta = summ["samples"][gsm]
    is_mt = gsm in {"GSM6481639", "GSM6481640", "GSM6481641"}
    # library-metric rows
    for ri, (col, ylab, _) in enumerate(rows):
        ax = axes[ri, ci]
        for f in ("top", "right"):
            ax.spines[f].set_visible(False)
        v = np.log10(np.maximum(sub[col], 1.0))
        vk = np.log10(np.maximum(kept[col], 1.0))
        ve = np.log10(np.maximum(excl[col], 1.0))
        if len(ve) > 0:
            ve = ve[~np.isnan(ve)]
        try:
            parts = ax.violinplot(vk, positions=[0], widths=0.85, showextrema=False)
            parts["bodies"][0].set_facecolor(KEPT_C)
            parts["bodies"][0].set_edgecolor("#5a5a5a")
            parts["bodies"][0].set_linewidth(0.5)
        except Exception:
            pass
        ax.scatter(rng.normal(0, 0.06, len(ve)), ve, s=1.2, color=EXCL_C,
                   alpha=0.45, linewidths=0, rasterized=True)
        ax.set_xlim(-1.4, 1.4)
        if ci == 0:
            ax.set_ylabel(ylab)
        else:
            ax.set_yticklabels([])
        ax.set_xticks([])
        ax.set_ylim(2.0, 4.7 if col == "n_genes" else 6.9)
        if ri == 0:
            ax.set_title(f"{labels[ci]}", loc="left", fontsize=9, weight="bold")
            ax.text(0.02, 0.97, f"kept {int(kept.shape[0]):,} / raw {int(sub.shape[0]):,}",
                    transform=ax.transAxes, va="top", ha="left",
                    fontsize=6.5, color="#555555")
    # filter-metric row: mixes two metrics (mito% for MT-row samples,
    # doublet score for MT-free). Each panel labels its own metric on the x
    # axis and keeps its own numeric ticks — a shared row-level y label would
    # misattribute the metric across columns.
    ax = axes[2, ci]
    for f in ("top", "right"):
        ax.spines[f].set_visible(False)
    if is_mt:
        v = sub["mt_frac"] * 100.0
        vk = kept["mt_frac"] * 100.0
        ve = excl["mt_frac"] * 100.0
        try:
            parts = ax.violinplot(vk, positions=[0], widths=0.85, showextrema=False)
            parts["bodies"][0].set_facecolor(KEPT_C)
            parts["bodies"][0].set_edgecolor("#5a5a5a")
            parts["bodies"][0].set_linewidth(0.5)
        except Exception:
            pass
        if len(ve):
            ax.scatter(rng.normal(0, 0.06, len(ve)), ve, s=1.4, color=EXCL_C,
                       alpha=0.5, linewidths=0, rasterized=True)
        ax.axhline(20.0, color=EXCL_C, lw=1.0, ls="--")
        ax.axhline(10.0, color="#444444", lw=0.8, ls=":")
        ax.set_ylim(0, 55)
        ax.text(0.02, 0.96,
                f"excl >20%: {int(meta['cells_mt20'])}",
                transform=ax.transAxes, va="top", ha="left", fontsize=6.5,
                color=EXCL_C)
        ax.text(0.02, 0.885,
                f"ref ≥10%: {int(meta['cells_mt10'])}",
                transform=ax.transAxes, va="top", ha="left", fontsize=6.5,
                color="#444444")
        ax.set_xlabel("mito %", fontsize=6.5)
    else:
        vs = sub["doublet_score"]
        vk = kept["doublet_score"]
        ve = excl["doublet_score"].dropna()
        try:
            parts = ax.violinplot(vk, positions=[0], widths=0.85, showextrema=False)
            parts["bodies"][0].set_facecolor(KEPT_C)
            parts["bodies"][0].set_edgecolor("#5a5a5a")
            parts["bodies"][0].set_linewidth(0.5)
        except Exception:
            pass
        if len(ve):
            ax.scatter(rng.normal(0, 0.06, len(ve)), ve, s=1.6, color=EXCL_C,
                       alpha=0.6, linewidths=0, rasterized=True)
        thr = meta["doublet_threshold"]
        if thr is not None:
            ax.axhline(thr, color=EXCL_C, lw=1.0, ls="--")
        ax.set_ylim(0, 1.0)
        n_dbl = int(meta["cells_doublet"])
        ax.text(0.02, 0.96, f"pred: {n_dbl}",
                transform=ax.transAxes, va="top", ha="left", fontsize=6.5,
                color=EXCL_C if n_dbl else "#444444")
        if thr is not None:
            ax.text(0.02, 0.885, f"thr: {thr:.2f}",
                    transform=ax.transAxes, va="top", ha="left", fontsize=6.5,
                    color="#444444")
        ax.set_xlabel("doublet score", fontsize=6.5)
    ax.set_xticks([])
    ax.set_xlim(-1.4, 1.4)

# shared legend (kept / excluded) hand-drawn in the reserved bottom band:
# deterministic coordinates beat fig.legend's box layout here.
from matplotlib.patches import Rectangle
fig.patches.append(Rectangle((0.432, 0.030), 0.014, 0.014,
                             facecolor=KEPT_C, edgecolor="#5a5a5a", linewidth=0.5,
                             transform=fig.transFigure))
fig.text(0.452, 0.034, "kept", fontsize=7, va="center", ha="left")
fig.patches.append(Rectangle((0.530, 0.030), 0.014, 0.014,
                             facecolor=EXCL_C, edgecolor=EXCL_C, linewidth=0.5,
                             transform=fig.transFigure))
fig.text(0.550, 0.034, "excluded", fontsize=7, va="center", ha="left")

# collision scan (figure-style prove-the-render, step 1)
rend = fig.canvas.get_renderer()
texts = [(t, t.get_window_extent(rend)) for t in fig.findobj(matplotlib.text.Text)
         if t.get_text().strip() and t.get_visible()]
frames = [(s, s.get_window_extent(rend)) for axx in fig.axes
          for s in axx.spines.values() if s.get_visible()]
own = {axx: set(axx.get_xticklabels(which="both") + axx.get_yticklabels(which="both"))
       for axx in fig.axes}
hits = [(a.get_text(), b.get_text()) for i, (a, ba) in enumerate(texts)
        for b, bb in texts[i + 1:] if ba.overlaps(bb)]
hits += [(t.get_text(), s) for t, bt in texts for s, bs in frames
         if bt.overlaps(bs) and t not in own[s.axes]]
assert not hits, f"collision scan failed: {hits}"

os.makedirs("figures", exist_ok=True)
fig.savefig(FIG, dpi=300, bbox_inches="tight")
print(f"saved {FIG}", flush=True)
