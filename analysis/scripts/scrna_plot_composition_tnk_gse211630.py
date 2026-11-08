# -*- coding: utf-8 -*-
"""Stage B figures:
Fig B1 geo_gse211630_composition: (a) family stacked composition per sample
   (b) subtype fraction shift (log10), sample dots + median ticks, star =
   sample-permutation p<=0.1 on the shift (min possible p given 5 samples).
Fig B2 geo_gse211630_tnk_function: T/NK subtype (rows) x functional signature
   (columns) matrix of per-sample mean module scores, NEOLP blue / EOLP orange;
   hollow circle = subtype has <10 cells in that sample; star = sample-perm
   p<=0.1 for that subtype x signature shift.

Reads the Stage-B tables; renders PNG (300 dpi) + PDF + 150-dpi QA copies and
runs the figure-style collision scan before success.

Panel API (2026-09-23)
----------------------
Refactored for figure recombination (Plan 3): the drawing bodies are now panel
functions that take a caller-supplied Axes, so a composing script can place the
same panel in a new plate. Geometry that is axes-relative (y-label offsets,
permutation stars, per-panel x limits) is unchanged; the figure-level suptitle
and footnote text stay with the original plates. `main()` still writes both
canonical figures, verified byte-identical to the pre-refactor output.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.text  # noqa: F401  (used by collision_scan)
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
os.chdir(ROOT)
TAB = os.path.join("results", "tables")
FIGDIR = "figures"
QA = os.path.join("analysis", "scrna", "_qa")

COND = ["NEOLP", "EOLP"]
COL = {"NEOLP": "#4477AA", "EOLP": "#EE7733"}
STAR = "#9A2A2A"
FAMILY_ORDER = ["T/NK", "B/Plasma", "Myeloid", "Mast", "Epithelial",
                "Fibroblast", "Pericyte/SMC", "Endothelial"]
FAM_COL = {"T/NK": "#882255", "B/Plasma": "#CC6677", "Myeloid": "#DDCC77",
           "Mast": "#117733", "Epithelial": "#88CCEE", "Fibroblast": "#44AA99",
           "Pericyte/SMC": "#999933", "Endothelial": "#332288"}
SHORT = {0: "c0 Fibroblast SFRP2+", 1: "c1 Fibroblast CXCL14+", 2: "c2 Pericyte/SMC",
         3: "c3 T resting", 4: "c4 T/NK stressed", 5: "c5 T/NK low-complex.",
         6: "c6 CD8+ GZMK+ EM", 7: "c7 CD4+ Treg", 8: "c8 T naive IL7R+",
         9: "c9 Plasmablast GZMB+", 10: "c10 Macrophage LYZ+",
         11: "c11 Epith. basal/prolif.", 13: "c13 Mast", 14: "c14 B cell",
         15: "c15 Plasma cell", 16: "c16 CD8+ GZMB+/CXCL13+",
         17: "c17 Endo lymph.", 19: "c19 Endo blood", 20: "c20 Myeloid S100-hi"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7,
                     "axes.linewidth": 0.6, "svg.fonttype": "none"})

comp_ps = pd.read_csv(os.path.join(TAB, "geo_gse211630_composition_per_sample.csv"))
comp_sh = pd.read_csv(os.path.join(TAB, "geo_gse211630_composition_shift.csv"))
mod_ps = pd.read_csv(os.path.join(TAB, "geo_gse211630_tnk_module_per_sample.csv"))
mod_sh = pd.read_csv(os.path.join(TAB, "geo_gse211630_tnk_module_shift.csv"))
samples = [f"NEOLP-{i}" for i in (1, 2, 3)] + [f"EOLP-{i}" for i in (1, 2)]
SIG_ORDER = list(dict.fromkeys(mod_ps.signature))       # preserves sorted order
cl_of_sub = dict(zip(mod_ps.subtype, mod_ps.cluster))
tnk_rows = sorted(mod_ps.subtype.unique(), key=lambda t: cl_of_sub[t])


def panel_tag(ax, tag, dx=-0.02, dy=1.02, fs=10.5):
    """Bold uppercase panel tag above the axes' top-left corner (matches the
    in-image (A)/(B) lettering used by main Figures 1-7).

    A falsy tag draws nothing, which lets a recomposition place the letter itself (so it can be
    sized and positioned by the plate's own hierarchy).
    """
    if not tag:
        return
    ax.text(dx, dy, tag, transform=ax.transAxes, fontsize=fs,
            fontweight="bold", ha="left", va="bottom")


def collision_scan(fig, rend):
    fig.canvas.draw()          # sync canvas renderer dpi with fig.dpi so window
                               # extents are measured consistently (savefig draws
                               # at other dpis would otherwise corrupt geometry)
    texts = [(t, t.get_window_extent(rend)) for t in fig.findobj(matplotlib.text.Text)
             if t.get_text().strip() and t.get_visible()]
    frames = [(s, s.get_window_extent(rend)) for ax in fig.axes
              for s in ax.spines.values() if s.get_visible()]
    own_ticks = {ax: set(ax.get_xticklabels(which="both") + ax.get_yticklabels(which="both"))
                 for ax in fig.axes}

    def hit(ba, bb):
        ox = min(ba.x1, bb.x1) - max(ba.x0, bb.x0)
        oy = min(ba.y1, bb.y1) - max(ba.y0, bb.y0)
        return ox > 2.5 and oy > 2.5     # ignore <2.5 px adjacency

    hits = [(a, b) for i, (a, ba) in enumerate(texts)
            for b, bb in texts[i + 1:] if hit(ba, bb)]
    hits += [(t, s) for t, bt in texts for s, bs in frames
             if hit(bt, bs) and t not in own_ticks[s.axes]]
    assert not hits, f"collision: {hits[:6]}"


# ======================================================================
# Fig B1 — panels
# ======================================================================
def panel_composition_bars(ax, tag=True, letter="A", cond_y=1.09):
    """(A) family stacked composition per sample."""
    fam_ps = comp_ps[comp_ps.level == "family"].copy()
    fam_ps["sample"] = pd.Categorical(fam_ps["sample"], categories=samples, ordered=True)
    fam_ps = fam_ps.sort_values("sample")
    bottom = np.zeros(len(samples))
    fam_shift_perm = dict(zip(comp_sh[comp_sh.level == "family"].group,
                              comp_sh[comp_sh.level == "family"].perm_p_sample))
    for f in FAMILY_ORDER:
        sub = fam_ps[fam_ps.group == f]
        if sub.empty:
            continue
        vals = sub.set_index("sample")["frac"].reindex(samples).fillna(0.0).values
        ax.bar(range(len(samples)), vals, bottom=bottom, width=0.62,
               color=FAM_COL[f], edgecolor="white", linewidth=0.5,
               label=(f + ("  *" if (fam_shift_perm.get(f, 1) or 1) <= 0.1 else "")))
        bottom += vals
    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, fontsize=6.8)
    ax.set_ylabel("fraction of annotated cells", fontsize=7.2, labelpad=14)
    ax.yaxis.set_label_coords(-0.18, 0.5)
    ax.set_ylim(0, 1.13)
    ax.set_title("Cell-type composition (families, per sample)", fontsize=8.4,
                 loc="left", pad=4)
    panel_tag(ax, letter, dx=-0.07, dy=1.0)
    totals = fam_ps.groupby("sample").n_cells.sum().reindex(samples).values
    for i, t in enumerate(totals):
        ax.text(i, 1.045, f"n={t:,}", ha="center", va="bottom", fontsize=5.6,
                color="#555555")
    # The y of the between-condition labels: a recomposition that raises the axis limit needs them
    # higher, otherwise they land on the per-sample n= labels drawn at y=1.045.
    ax.text(1.0, cond_y, "NEOLP (n=3)", ha="center", fontsize=6.6, color=COL["NEOLP"])
    ax.text(3.5, cond_y, "EOLP (n=2)", ha="center", fontsize=6.6, color=COL["EOLP"])
    ax.axvline(2.5, color="#BBBBBB", lw=0.8, ls="--")
    # family legend below panel (a) in the empty bottom margin (2 rows x 4 cols);
    # a tall outside legend between the two panels would collide with panel (b)'s
    # top row labels.
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.11), ncol=4,
              frameon=False, fontsize=6.0, columnspacing=1.4, handletextpad=0.5,
              labelspacing=0.6)
    return ax


def panel_subtype_dots(ax, tag=True, letter="B"):
    """(B) subtype fraction dots (log10), sample dots + median ticks + perm stars."""
    sub_ps = comp_ps[comp_ps.level == "subtype"].copy()
    sub_ps["cluster"] = sub_ps["cluster"].astype(int)
    # subtype short label keyed by cluster; family keyed by cluster id.
    # NOTE: at the subtype level `group` holds the descriptive subtype name
    # (e.g. "CD8+ T (GZMK+ effector-memory)"), not the family -- the cluster ->
    # family mapping is taken from the shift table, which carries both columns.
    sub_sh = comp_sh[comp_sh.level == "subtype"]
    cl2fam = {int(r.cluster): r.family for _, r in sub_sh.iterrows()}
    fam_of_cl = {c: cl2fam[c] for c in sub_ps.cluster.unique()}
    lab_of_cl = {c: SHORT[c] for c in sub_ps.cluster.unique()}
    minpos = sub_ps[sub_ps.frac > 0].frac.min()
    floor = minpos / 2
    sub_sh_perm = dict(zip(sub_sh.group, sub_sh.perm_p_sample))

    ypos = {}
    y = 0
    ytick, yticklab, ytickcol, seps = [], [], [], []
    for f in FAMILY_ORDER:
        cs = sorted(c for c, fam in fam_of_cl.items() if fam == f)
        for c in cs:
            g = sub_ps[sub_ps.cluster == c].group.iloc[0]
            ypos[g] = y
            ytick.append(y)
            yticklab.append(lab_of_cl[c])
            rgb = np.asarray(matplotlib.colors.to_rgb(FAM_COL[f]))
            ytickcol.append(matplotlib.colors.to_hex(np.clip(rgb * 0.55 + 0.05, 0, 1)))
            y += 1
        seps.append(y - 0.5)
    ylim_top = y - 1

    for g, yy in ypos.items():
        row = sub_ps[sub_ps.group == g]
        for cn_ in COND:
            for _, q in row[row.condition == cn_].iterrows():
                x = np.log10(q.frac) if q.frac > 0 else np.log10(floor)
                ax.scatter(x, yy, s=22, color=COL[cn_], edgecolor="white", lw=0.4,
                           zorder=3, marker="o" if q.frac > 0 else "s")
            v = row[row.condition == cn_].frac.values
            med = np.log10(np.median(v[v > 0])) if (v > 0).any() else np.log10(floor)
            ax.plot([med, med], [yy - 0.28, yy + 0.28], color=COL[cn_], lw=1.6, zorder=4)
    ax.set_yticks(ytick)
    ax.set_yticklabels(yticklab, fontsize=6.2)
    for lab, colr in zip(ax.get_yticklabels(), ytickcol):
        lab.set_color(colr)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.8, ylim_top + 0.8)
    ax.set_xlim(np.log10(floor) - 0.55, np.log10(sub_ps.frac.max()) + 0.35)
    xmin, xmax = ax.get_xlim()
    sx = xmin + 0.05 * (xmax - xmin)
    for g, yy in ypos.items():
        if sub_sh_perm.get(g, 1) <= 0.1:
            ax.text(sx, yy, "*", color=STAR, fontsize=8.5, ha="left", va="center", zorder=5)
    for sp in seps:
        ax.axhline(sp, color="#CCCCCC", lw=0.6, zorder=1)
    ax.set_xlabel(r"fraction of annotated cells in the sample  (log$_{10}$ scale)",
                  fontsize=7.0, labelpad=10)
    xt = np.arange(-4, 0, 1)
    ax.set_xticks(xt)
    ax.set_xticklabels([f"$10^{{{int(t)}}}$" for t in xt], fontsize=6.2)
    ax.set_title("Subtype fractions in NEOLP vs EOLP", fontsize=8.4, loc="left", pad=4)
    panel_tag(ax, letter, dx=-0.07, dy=1.0)
    return ax


# ======================================================================
# Fig B2 — panel (one signature column)
# ======================================================================
def panel_tnk_signature(ax, sig, ylabels=False):
    """One T/NK module-score panel for `sig` (rows = T/NK subtypes)."""
    for row_i, t in enumerate(tnk_rows):
        # row order matches the y-tick labels (ascending cluster id:
        # c3 at tick 0/bottom -> c16 at top).  Previously yy was
        # len-1-row_i, which inverted the data rows relative to the
        # labels (geometry-verified); fixed to yy = row_i.
        yy = row_i
        sub = mod_ps[(mod_ps.signature == sig) & (mod_ps.subtype == t)]
        for cn_ in COND:
            part = sub[sub.condition == cn_]
            for _, q in part.iterrows():
                if q.n_cells >= 10:
                    ax.scatter(q.mean_score, yy, s=30, color=COL[cn_],
                               edgecolor="white", linewidth=0.4, zorder=3)
                else:                                   # sparse sample (<10 cells)
                    ax.scatter(q.mean_score, yy, s=28, facecolor="none",
                               edgecolor=COL[cn_], linewidth=1.2, zorder=3)
            if len(part):
                med = part.mean_score.median()
                # median tick under the sample dots (zorder 2 < 4): group
                # median must never occlude a sample marker lying at the median
                ax.plot([med, med], [yy - 0.30, yy + 0.30], color=COL[cn_],
                        lw=1.6, zorder=2)
    ax.margins(x=0.18, y=0.12)
    # prune='both': MaxNLocator's edge ticks may fall outside the data limits
    # (e.g. tick 2.4 vs xlim 2.0); unclipped labels would leak into the
    # inter-panel gaps and collide with the neighbouring panel's labels.
    ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=4, prune="both"))
    ax.set_yticks(range(len(tnk_rows)))
    ax.set_yticklabels([])
    ax.set_ylim(-0.6, len(tnk_rows) - 0.4)
    ax.set_title(sig, fontsize=7.6, pad=6)
    ax.tick_params(axis="y", length=0, labelsize=5.8)
    ax.tick_params(axis="x", labelsize=5.8)
    ax.set_xlabel("module score", fontsize=6.4, labelpad=5)
    # sample-level permutation star: placed in the upper whitespace of the row
    # band, vertically above the highest point of the higher-median group --
    # never on a marker (rows are ~70 px tall; dots span only the central
    # ~20 px), never in the inter-panel gutter where it could be misattributed
    # to the neighbouring panel.
    xspan = ax.get_xlim()[1] - ax.get_xlim()[0]
    for row_i, t in enumerate(tnk_rows):
        msh = mod_sh[(mod_sh.signature == sig) & (mod_sh.subtype == t)]
        if len(msh) and msh.perm_p_sample.iloc[0] <= 0.1:
            cell = mod_ps[(mod_ps.signature == sig) & (mod_ps.subtype == t)]
            meds = {cn_: cell.loc[cell.condition == cn_, "mean_score"].median()
                    for cn_ in COND}
            meds = {k: v for k, v in meds.items() if not pd.isna(v)}
            if meds:
                hi = max(meds, key=meds.get)
                hi_max = cell.loc[cell.condition == hi, "mean_score"].max()
                ax.text(hi_max + 0.01 * xspan, row_i + 0.34,
                        "*", color=STAR, fontsize=8.0, ha="center", va="center",
                        zorder=5)
    for s_ in ax.spines.values():
        s_.set_linewidth(0.6)
    if ylabels:
        ax.set_yticklabels([SHORT[int(cl_of_sub[t])] for t in tnk_rows], fontsize=6.4)
        ax.tick_params(axis="y", length=0)
    return ax


TNK_LEGEND = [("NEOLP sample", dict(marker="o", ls="", color=COL["NEOLP"], ms=5.4)),
              ("EOLP sample", dict(marker="o", ls="", color=COL["EOLP"], ms=5.4)),
              ("<10 cells of subtype",
               dict(marker="o", ls="", mfc="none", mec="#666666", ms=4.8))]


def tnk_legend_handles():
    return [Line2D([], [], label=lab, **kw) for lab, kw in TNK_LEGEND]


# ======================================================================
# Original plates
# ======================================================================
def build_b1():
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 6.6),
                             gridspec_kw={"width_ratios": [0.9, 1.9],
                                          "wspace": 0.42, "left": 0.06, "right": 0.995,
                                          "top": 0.90, "bottom": 0.20})
    axa, axb = axes
    panel_composition_bars(axa)
    panel_subtype_dots(axb)
    fig.text(0.995, 0.02, "* sample-level permutation p≤0.1 (min possible with 5 samples, two-sided\n"
             "over the 10 label splits).  Blue = NEOLP (3 samples), orange = EOLP (2);\n"
             "square = fraction 0 (drawn at the floor).  Row colour = family (as stacked left).",
             fontsize=5.8, color="#333333", ha="right", va="bottom")
    fig.suptitle("OLP composition: non-erosive vs erosive samples (GSE211630, n=46,446 annotated cells)",
                 fontsize=9.2, weight="bold", x=0.28, ha="left")
    return fig


def build_b2():
    fig, axs = plt.subplots(1, len(SIG_ORDER), figsize=(12.2, 4.9), sharey=False,
                            gridspec_kw={"left": 0.115, "right": 0.985, "top": 0.84,
                                         "bottom": 0.175, "wspace": 0.42})
    for i, (ax, sig) in enumerate(zip(axs, SIG_ORDER)):
        panel_tnk_signature(ax, sig, ylabels=(i == 0))
    fig.legend(handles=tnk_legend_handles(), loc="lower center", ncol=3, frameon=False,
               fontsize=6.4, bbox_to_anchor=(0.30, 0.012), handletextpad=0.6,
               columnspacing=1.4)
    fig.text(0.985, 0.018, "* sample-level permutation p≤0.1 for the NEOLP−EOLP shift "
             "(two-sided, 10 label splits).  Vertical tick = group median of the "
             "per-sample module scores;\nhollow marker = <10 cells of the subtype in "
             "that sample (marker colour still encodes the group).  x-axis scale is "
             "per panel.", fontsize=6.0, color="#333333", ha="right", va="bottom")
    fig.suptitle("T/NK functional states (module scores of marker signatures), NEOLP vs EOLP",
                 fontsize=9.4, weight="bold", x=0.24, ha="left")
    return fig


def main():
    os.makedirs(QA, exist_ok=True)
    fig = build_b1()
    fig.savefig(os.path.join(FIGDIR, "geo_gse211630_composition.png"), dpi=300)
    fig.savefig(os.path.join(FIGDIR, "geo_gse211630_composition.pdf"))
    fig.savefig(os.path.join(QA, "fig_composition.png"), dpi=150)
    collision_scan(fig, fig.canvas.get_renderer())
    plt.close(fig)
    print("Fig B1 written: figures/geo_gse211630_composition.png")

    fig = build_b2()
    fig.savefig(os.path.join(FIGDIR, "geo_gse211630_tnk_function.png"), dpi=300)
    fig.savefig(os.path.join(FIGDIR, "geo_gse211630_tnk_function.pdf"))
    fig.savefig(os.path.join(QA, "fig_tnk_function.png"), dpi=150)
    collision_scan(fig, fig.canvas.get_renderer())
    plt.close(fig)
    print("Fig B2 written: figures/geo_gse211630_tnk_function.png")


if __name__ == "__main__":
    main()
