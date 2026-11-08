# -*- coding: utf-8 -*-
"""Supplementary Figure 7: lineage-marker dot plot used for GSE211630 atlas annotation.

gse211630_olp_lognorm.h5ad was built with ad.concat(join="inner") and is
17,669 genes; 17 of the curated markers are absent from that intersection.
This plate therefore scores every marker on the UNION of the five OLP samples
(data/processed/geo_GSE211630/ann/GSM648163{8..42}_filtered.h5ad) while taking
the cluster labels from the merged object. Markers absent from all five
matrices are dropped (TRAC); a coverage strip below the dot plot gives, per
marker, the number of the five samples whose matrix contains the gene.

Outputs
-------
figures/supplementary/supp_fig07_marker_dotplot.png (300 dpi) + .pdf

Panel API (2026-09-23)
----------------------
Refactored for figure recombination (Plan 3): the drawing bodies are now panel functions taking
caller-supplied Axes. The five-matrix aggregation is cached by `load()` so a composing script pays
for it once. Drawing statements are unchanged; `main()` still writes the canonical figure, verified
byte-identical to the pre-refactor output.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

exec(compile(open(r"\\?\D:\wisp-science\skills\figure-style\runtime.py",
                  encoding="utf-8").read(), "runtime.py", "exec"))

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
os.chdir(ROOT)
OUT = os.path.join("data", "processed", "geo_GSE211630")
ANN = os.path.join(OUT, "ann")
SCAN = os.path.join(OUT, "scanpy")
FIGDIR = os.path.join("figures", "supplementary")
os.makedirs(FIGDIR, exist_ok=True)

apply_figure_style(sizes=(7, 6.3, 5.7))

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
OLP = [("GSM6481638", "NEOLP-1"), ("GSM6481639", "NEOLP-2"),
       ("GSM6481640", "NEOLP-3"), ("GSM6481641", "EOLP-1"),
       ("GSM6481642", "EOLP-2")]
LEID = "leiden_0.5"

import anndata as ad
import scanpy as sc

_CACHE = None


def load():
    """Aggregate marker expression per cluster across the five OLP matrices (cached).

    Runs once per process: it reads five h5ad files, so a composing script that uses this panel
    should not pay for it twice.
    """
    global _CACHE
    if _CACHE is not None:
        return _CACHE

    olp = ad.read_h5ad(os.path.join(SCAN, "gse211630_olp_lognorm.h5ad"))
    leiden = olp.obs[LEID].astype(int).values
    slab = olp.obs["label"].astype(str).values
    clu = sorted(np.unique(leiden))

    gene_all = [g for genes in MARKERS.values() for g in genes]
    gi = {g: k for k, g in enumerate(gene_all)}
    n_g = len(gene_all)
    acc_n = np.zeros((len(clu), n_g))
    acc_sum = np.zeros((len(clu), n_g))
    acc_pos = np.zeros((len(clu), n_g))
    cov = np.zeros(n_g, dtype=int)

    for gsm, lab in OLP:
        a = sc.read_h5ad(os.path.join(ANN, f"{gsm}_filtered.h5ad"))
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
        block = leiden[slab == lab]
        if len(block) != a.n_obs:
            raise RuntimeError(f"{lab}: {a.n_obs} cells vs {len(block)} cluster labels")
        have = [g for g in gene_all if g in a.var_names]
        Xs = a[:, have].X
        Xs = np.asarray(Xs.todense()) if hasattr(Xs, "todense") else np.asarray(Xs)
        for j, g in enumerate(have):
            k = gi[g]
            cov[k] += 1
            col = Xs[:, j]
            for ci, c in enumerate(clu):
                m = block == c
                n = int(m.sum())
                if n == 0:
                    continue
                v = col[m]
                acc_n[ci, k] += n
                acc_sum[ci, k] += float(v.sum())
                acc_pos[ci, k] += float((v > 0).sum())
        del a, Xs

    raw_mean = np.divide(acc_sum, acc_n, out=np.zeros_like(acc_sum), where=acc_n > 0)
    raw_frac = np.divide(acc_pos, acc_n, out=np.zeros_like(acc_pos), where=acc_n > 0)

    dropped = [g for g in gene_all if cov[gi[g]] == 0]
    print("markers absent from ALL five OLP matrices (dropped):", dropped, flush=True)

    keep = [g for g in gene_all if cov[gi[g]] > 0]
    kidx = [gi[g] for g in keep]
    gene_order = keep
    mean_exp = raw_mean[:, kidx]
    frac_exp = raw_frac[:, kidx]
    cov_p = cov[kidx]
    print(f"markers plotted: {len(gene_order)}/{len(gene_all)} "
          f"| coverage min {cov_p.min()} max {cov_p.max()} of 5", flush=True)
    print("coverage per marker:", {g: int(c) for g, c in zip(gene_order, cov_p)},
          flush=True)

    mu = mean_exp.mean(axis=0)
    sd = mean_exp.std(axis=0)
    z = (mean_exp - mu) / np.where(sd > 0, sd, 1.0)
    print("NaN cells in z:", int(np.isnan(z).sum()), "of", z.size, flush=True)

    _CACHE = dict(clu=clu, gene_order=gene_order, z=z, frac_exp=frac_exp, cov_p=cov_p,
                  n_g=n_g, dropped=dropped)
    return _CACHE


def family_bounds(gene_order):
    """[(family, first index, last index)] over the plotted marker order."""
    bounds, k = [], 0
    for g, genes in MARKERS.items():
        n = sum(1 for x in genes if x in set(gene_order))
        if n:
            bounds.append((g, k, k + n - 1))
        k += n
    return bounds


def panel_dotplot(ax):
    """Marker x cluster dot plot. Returns the scatter artist (for a caller-placed colorbar)."""
    D = load()
    clu, gene_order, z, frac_exp = D["clu"], D["gene_order"], D["z"], D["frac_exp"]

    bounds = family_bounds(gene_order)
    pos = np.arange(len(gene_order))
    for _g, lo, _hi in bounds:
        ax.axvline(lo - 0.5, color="#E0E0E0", lw=0.5)

    vmax = max(1.0, float(np.abs(np.nan_to_num(z)).max()) * 0.8)
    ok = ~np.isnan(z)
    scat = ax.scatter(np.repeat(pos[None, :], len(clu), axis=0),
                      np.repeat(np.array(clu)[:, None], len(gene_order), axis=1),
                      c=np.nan_to_num(z), cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                      s=np.where(ok, 2.0 + 78.0 * np.clip(frac_exp, 0, 1) ** 0.5, 0.0),
                      linewidths=0, alpha=0.9, rasterized=True)
    ax.set_xticks(pos)
    ax.set_xticklabels(gene_order, rotation=90, fontsize=5.5)
    ax.set_yticks(clu)
    ax.set_ylim(clu[0] - 0.6, clu[-1] + 2.7)
    ax.set_xlim(-0.5, len(gene_order) - 0.5)
    ax.set_ylabel("Leiden cluster")
    set_frame(ax)

    for i, (g, lo, hi) in enumerate(bounds):
        ax.text((lo + hi) / 2.0, len(clu) - 0.45 + (1.35 if i % 2 else 0.0), g,
                ha="center", va="bottom", fontsize=5.2, color="#444444")
    return scat


def panel_coverage(axc):
    """Coverage strip: number of the five OLP matrices containing each marker."""
    D = load()
    gene_order, cov_p = D["gene_order"], D["cov_p"]
    pos = np.arange(len(gene_order))
    axc.bar(pos, cov_p, width=0.72, color="#B0B0B0", edgecolor="none")
    axc.axhline(5, color="#8C8C8C", lw=0.5, ls=":")
    axc.set_xlim(-0.5, len(gene_order) - 0.5)
    axc.set_ylim(0, 5.4)
    axc.set_yticks([1, 3, 5])
    axc.set_ylabel("# samples", fontsize=5.5, labelpad=1.5)
    axc.set_xticks([])
    axc.tick_params(axis="y", labelsize=5.2, length=2, pad=1.5)
    set_frame(axc)
    return axc


def add_colorbar(fig, scat, rect):
    """Horizontal z-scale colorbar on its own axes."""
    cax = fig.add_axes(rect)
    cbar = fig.colorbar(scat, cax=cax, orientation="horizontal")
    cbar.set_label("z(mean log1p expr.)", fontsize=5.7, labelpad=1.5)
    cbar.ax.tick_params(labelsize=5.2, length=2, pad=1.5)
    return cbar


def frac_legend_handles():
    """Proxy handles encoding the dot area = fraction of expressing cells."""
    prox = []
    for f in (0.05, 0.25, 0.5, 0.75, 1.0):
        s = 2.0 + 78.0 * np.sqrt(f)
        prox.append(plt.Line2D([], [], marker="o", ls="", color="#8C8C8C",
                               markersize=2 * np.sqrt(s / np.pi), label=f"{f:g}"))
    return prox


def add_frac_legend(fig, at=(0.525, 0.088)):
    fig.legend(handles=frac_legend_handles(), loc="center left", bbox_to_anchor=at,
               frameon=False, ncol=5, title="frac. expressing", title_fontsize=5.7,
               fontsize=5.7, columnspacing=1.2, handletextpad=0.35,
               borderpad=0.0, labelspacing=0.5)


def collision_scan(fig, ax, axc, cbar):
    """Text-collision / off-canvas QA, plus the gene-label vs coverage-strip check."""
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    arts = [t for t in ax.texts if t.get_text().strip()]
    for a_ in (ax, axc):
        for t in (a_.xaxis.label, a_.yaxis.label):
            if t.get_text().strip():
                arts.append(t)
        if a_.get_xticks().size:
            arts += [t for t in a_.get_xticklabels() if t.get_text().strip()]
        if a_.get_yticks().size:
            arts += [t for t in a_.get_yticklabels() if t.get_text().strip()]
    if cbar.ax.xaxis.label.get_text().strip():
        arts.append(cbar.ax.xaxis.label)
    arts += [t for t in cbar.ax.get_xticklabels() if t.get_text().strip()]
    for leg in fig.legends:
        arts += [t for t in leg.get_texts() if t.get_text().strip()]
        if leg.get_title().get_text().strip():
            arts.append(leg.get_title())
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

    # the rotated gene labels must stay clear of the coverage strip beneath them:
    # text/text scanning cannot catch a label running into a neighbouring axes.
    bxc = axc.get_window_extent(renderer=rend)
    intr = []
    for t in ax.get_xticklabels():
        if not t.get_text().strip():
            continue
        b = t.get_window_extent(renderer=rend)
        if (min(b.x1, bxc.x1) - max(b.x0, bxc.x0) > TOL and
                min(b.y1, bxc.y1) - max(b.y0, bxc.y0) > TOL):
            intr.append(t.get_text())
    print("gene labels intruding into the coverage strip:", intr, flush=True)
    assert not intr, "collision scan failed (gene labels vs coverage strip)"


def build():
    """Original 180 x 112 mm plate: dot plot, coverage strip, colorbar and dot-area legend."""
    MM = 1 / 25.4
    fig = plt.figure(figsize=(180 * MM, 112 * MM))
    # the rotated gene labels need ~15 mm of clear space, so the coverage strip sits
    # well below the main axes instead of directly beneath it
    ax = fig.add_axes([0.050, 0.410, 0.940, 0.500])
    axc = fig.add_axes([0.050, 0.195, 0.940, 0.043])
    scat = panel_dotplot(ax)
    panel_coverage(axc)
    cbar = add_colorbar(fig, scat, [0.108, 0.075, 0.30, 0.027])
    add_frac_legend(fig)
    return fig, ax, axc, cbar


def main():
    fig, ax, axc, cbar = build()
    collision_scan(fig, ax, axc, cbar)
    fig.savefig(os.path.join(FIGDIR, "supp_fig07_marker_dotplot.png"), dpi=300,
                bbox_inches="tight")
    fig.savefig(os.path.join(FIGDIR, "supp_fig07_marker_dotplot.pdf"),
                bbox_inches="tight")
    D = load()
    print("clusters:", len(D["clu"]), "| families:", len(family_bounds(D["gene_order"])),
          flush=True)
    print("saved supp_fig07_marker_dotplot.{png,pdf}", flush=True)


if __name__ == "__main__":
    main()
