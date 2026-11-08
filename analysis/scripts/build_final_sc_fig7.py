# -*- coding: utf-8 -*-
"""Supplementary figure = original Fig9 (composition) + Fig10 (T/NK module scores).

Moved out of the main figures because it is not part of the main MR storyline.
Output: figures/supplementary/supp_fig09_composition_tnk.{png,pdf}

Renumbered 2026-09-13: the supplementary set was reordered to follow the
manuscript narrative, so this plate moved from Supplementary Figure 6 to
Supplementary Figure 9 (three single-cell plates now precede it: the dataset/QC
plate, the atlas clustering plate and the marker dot plot). The plate content is
unchanged and the script filename is retained as the historical stage name.

Layout (12.2 x 11.4 in, 300 dpi PNG + vector PDF):
  top row:   A family stacked composition per sample  |  B subtype fraction dots
  bottom row: C T/NK functional module scores (7 signatures x 7 subtypes)
Panels re-tagged A/B/C. Reads the same Stage-B tables as
scrna_plot_composition_tnk_gse211630.py and reuses its drawing code verbatim.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
os.chdir(ROOT)
TAB = os.path.join("results", "tables")
# 2026-09: this plate is not part of the main-analysis storyline -> supplementary
FIGDIR = os.path.join("figures", "supplementary")
os.makedirs(FIGDIR, exist_ok=True)
QA = os.path.join("analysis", "scrna", "_qa")
os.makedirs(QA, exist_ok=True)

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


def panel_tag(ax, tag, dx=-0.02, dy=1.02, fs=10.5):
    ax.text(dx, dy, tag, transform=ax.transAxes, fontsize=fs,
            fontweight="bold", ha="left", va="bottom")


def collision_scan(fig, rend):
    fig.canvas.draw()
    texts = [(t, t.get_window_extent(rend)) for t in fig.findobj(matplotlib.text.Text)
             if t.get_text().strip() and t.get_visible()]
    frames = [(s, s.get_window_extent(rend)) for ax in fig.axes
              for s in ax.spines.values() if s.get_visible()]
    own_ticks = {ax: set(ax.get_xticklabels(which="both") + ax.get_yticklabels(which="both"))
                 for ax in fig.axes}

    def hit(ba, bb):
        ox = min(ba.x1, bb.x1) - max(ba.x0, bb.x0)
        oy = min(ba.y1, bb.y1) - max(ba.y0, bb.y0)
        return ox > 2.5 and oy > 2.5

    hits = [(a, b) for i, (a, ba) in enumerate(texts)
            for b, bb in texts[i + 1:] if hit(ba, bb)]
    hits += [(t, s) for t, bt in texts for s, bs in frames
             if hit(bt, bs) and t not in own_ticks[s.axes]]
    assert not hits, f"collision: {hits[:6]}"


# ------------------------------------------------------------------ figure
fig = plt.figure(figsize=(12.2, 11.4))
gs = fig.add_gridspec(2, 1, height_ratios=[1.32, 1.0],
                      left=0.115, right=0.995, top=0.985, bottom=0.135,
                      hspace=0.34)
gs_top = gs[0].subgridspec(1, 2, width_ratios=[0.9, 1.9], wspace=0.42)
axa = fig.add_subplot(gs_top[0])
axb = fig.add_subplot(gs_top[1])
gs_bot = gs[1].subgridspec(1, 7, wspace=0.42)
axs = [fig.add_subplot(gs_bot[i]) for i in range(7)]

# =====================================================================
# Panel A: family stacked composition per sample
# =====================================================================
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
    axa.bar(range(len(samples)), vals, bottom=bottom, width=0.62,
            color=FAM_COL[f], edgecolor="white", linewidth=0.5,
            label=(f + ("  *" if (fam_shift_perm.get(f, 1) or 1) <= 0.1 else "")))
    bottom += vals
axa.set_xticks(range(len(samples)))
axa.set_xticklabels(samples, fontsize=6.8)
axa.set_ylabel("fraction of annotated cells", fontsize=7.2, labelpad=14)
axa.yaxis.set_label_coords(-0.18, 0.5)
axa.set_ylim(0, 1.13)
axa.set_yticks(np.arange(0, 1.01, 0.2))
axa.set_title("Cell-type composition (families, per sample)", fontsize=8.4, loc="left", pad=4)
panel_tag(axa, "A", dx=-0.07, dy=1.0)
totals = fam_ps.groupby("sample").n_cells.sum().reindex(samples).values
for i, t in enumerate(totals):
    axa.text(i, 1.045, f"n={t:,}", ha="center", va="bottom", fontsize=5.6, color="#555555")
axa.text(1.0, 1.09, "NEOLP (n=3)", ha="center", fontsize=6.6, color=COL["NEOLP"])
axa.text(3.5, 1.09, "EOLP (n=2)", ha="center", fontsize=6.6, color=COL["EOLP"])
axa.axvline(2.5, color="#BBBBBB", lw=0.8, ls="--")
axa.legend(loc="upper center", bbox_to_anchor=(0.5, -0.11), ncol=4,
           frameon=False, fontsize=6.0, columnspacing=1.4, handletextpad=0.5,
           labelspacing=0.6)

# =====================================================================
# Panel B: subtype fraction dots
# =====================================================================
sub_ps = comp_ps[comp_ps.level == "subtype"].copy()
sub_ps["cluster"] = sub_ps["cluster"].astype(int)
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
            axb.scatter(x, yy, s=22, color=COL[cn_], edgecolor="white", lw=0.4,
                        zorder=3, marker="o" if q.frac > 0 else "s")
        v = row[row.condition == cn_].frac.values
        med = np.log10(np.median(v[v > 0])) if (v > 0).any() else np.log10(floor)
        axb.plot([med, med], [yy - 0.28, yy + 0.28], color=COL[cn_], lw=1.6, zorder=4)
axb.set_yticks(ytick)
axb.set_yticklabels(yticklab, fontsize=6.2)
for lab, colr in zip(axb.get_yticklabels(), ytickcol):
    lab.set_color(colr)
axb.tick_params(axis="y", length=0)
axb.set_ylim(-0.8, ylim_top + 0.8)
axb.set_xlim(np.log10(floor) - 0.55, np.log10(sub_ps.frac.max()) + 0.35)
xmin, xmax = axb.get_xlim()
sx = xmin + 0.05 * (xmax - xmin)
for g, yy in ypos.items():
    if sub_sh_perm.get(g, 1) <= 0.1:
        axb.text(sx, yy, "*", color=STAR, fontsize=8.5, ha="left", va="center", zorder=5)
for sp in seps:
    axb.axhline(sp, color="#CCCCCC", lw=0.6, zorder=1)
axb.set_xlabel(r"fraction of annotated cells in the sample  (log$_{10}$ scale)",
               fontsize=7.0, labelpad=10)
xt = np.arange(-4, 0, 1)
axb.set_xticks(xt)
axb.set_xticklabels([f"$10^{{{int(t)}}}$" for t in xt], fontsize=6.2)
axb.set_title("Subtype fractions in NEOLP vs EOLP", fontsize=8.4, loc="left", pad=4)
panel_tag(axb, "B", dx=-0.07, dy=1.0)

# =====================================================================
# Panel C: T/NK functional module scores (7 signatures)
# =====================================================================
SIG_ORDER = list(dict.fromkeys(mod_ps.signature))
cl_of_sub = dict(zip(mod_ps.subtype, mod_ps.cluster))
tnk_rows = sorted(mod_ps.subtype.unique(), key=lambda t: cl_of_sub[t])
for ax, sig in zip(axs, SIG_ORDER):
    for row_i, t in enumerate(tnk_rows):
        sub = mod_ps[(mod_ps.signature == sig) & (mod_ps.subtype == t)]
        for cn_ in COND:
            part = sub[sub.condition == cn_]
            for _, q in part.iterrows():
                if q.n_cells >= 10:
                    ax.scatter(q.mean_score, row_i, s=30, color=COL[cn_],
                               edgecolor="white", linewidth=0.4, zorder=3)
                else:
                    ax.scatter(q.mean_score, row_i, s=28, facecolor="none",
                               edgecolor=COL[cn_], linewidth=1.2, zorder=3)
            if len(part):
                med = part.mean_score.median()
                ax.plot([med, med], [row_i - 0.30, row_i + 0.30], color=COL[cn_], lw=1.6, zorder=2)
    ax.margins(x=0.18, y=0.12)
    ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=4, prune="both"))
    ax.set_yticks(range(len(tnk_rows)))
    ax.set_yticklabels([])
    ax.set_ylim(-0.6, len(tnk_rows) - 0.4)
    ax.set_title(sig, fontsize=7.6, pad=6)
    ax.tick_params(axis="y", length=0, labelsize=5.8)
    ax.tick_params(axis="x", labelsize=5.8)
    ax.set_xlabel("module score", fontsize=6.4, labelpad=5)
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
                ax.text(hi_max + 0.01 * xspan, row_i + 0.34, "*", color=STAR,
                        fontsize=8.0, ha="center", va="center", zorder=5)
    for s_ in ax.spines.values():
        s_.set_linewidth(0.6)
axs[0].set_yticklabels([SHORT[int(cl_of_sub[t])] for t in tnk_rows], fontsize=6.4)
axs[0].tick_params(axis="y", length=0)
panel_tag(axs[0], "C", dx=-0.055, dy=1.06)

# ------------------------------------------------------------------ legend + notes
h = [Line2D([], [], marker="o", ls="", color=COL["NEOLP"], ms=5.4, label="NEOLP sample"),
     Line2D([], [], marker="o", ls="", color=COL["EOLP"], ms=5.4, label="EOLP sample"),
     Line2D([], [], marker="o", ls="", mfc="none", mec="#666666", ms=4.8, label="<10 cells of subtype")]
fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, fontsize=6.4,
           bbox_to_anchor=(0.44, 0.062), handletextpad=0.6, columnspacing=1.4)

fig.text(0.995, 0.024,
         "* sample-level permutation p\u22640.1 (two-sided, 10 label splits).  "
         "A\u2013B: square = fraction 0 (floor); row colour = family.  "
         "C: vertical tick = group median.",
         fontsize=5.9, color="#333333", ha="right", va="bottom")
# (no figure-level title: this is a panel of a multi-figure plate, the caption
#  carries the heading; grid top raised from 0.955 to 0.985 to reclaim the band)

fig.savefig(os.path.join(FIGDIR, "supp_fig09_composition_tnk.png"), dpi=300)
fig.savefig(os.path.join(FIGDIR, "supp_fig09_composition_tnk.pdf"))
fig.savefig(os.path.join(QA, "fig7_merged.png"), dpi=150)
collision_scan(fig, fig.canvas.get_renderer())
plt.close(fig)
print("Supplementary fig (composition + T/NK) written: "
      "figures/supplementary/supp_fig09_composition_tnk.png")
