# -*- coding: utf-8 -*-
"""
Fig C1 (MR drug-target reverse mapping, GSE211630): EOLP-vs-NEOLP direction
of the five finalized OLP MR targets inside each of the 19 cell subtypes
(panel a) and the tissue-level decomposition of each target's expression
change into composition vs within-celltype vs interaction (panel b, horizontal
waterfall). Observed tissue polarity is read against the MR-implied polarity
(risk gene  => EOLP message expected up; protective gene => expected down).

Colour bindings (house palette):
  NEOLP #4477AA / EOLP #EE7733 (not shown here); family colours only tint the
  y tick labels in panel a (same palette as Fig B1).  Panel a marker colours:
  crimson = EOLP up (supported), navy = EOLP down (supported), light grey =
  no supported change, x = not detected; black ring = supported by BOTH the
  cell-level BH-FDR and the sample-permutation track.  Panel b segments:
  ochre = composition, steel = within-celltype, grey = interaction; black
  diamond = observed tissue delta; gene label carries the MR badge.

Outputs: figures/geo_gse211630_mr_backmap.png (+pdf, 300 dpi), QA copy 150 dpi.
Inputs : results/tables/geo_gse211630_mr_backmap_{gene,subtype}.csv
Run    : python analysis/scripts/plot_mr_backmap_gse211630.py  (matplotlib only)
"""

import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = os.getcwd()
TAB = os.path.join(ROOT, "results", "tables")
FIGDIR = os.path.join(ROOT, "figures")
QA = os.path.join(ROOT, "analysis", "scrna", "_qa")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(QA, exist_ok=True)

# ---------- figure style (skill runtime first, inline fallback otherwise)
try:
    apply_figure_style  # noqa: F821
    apply_figure_style(sizes=(8.2, 7.0, 6.2), frame="open")
except NameError:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8.2,
        "axes.linewidth": 0.6, "axes.spines.top": False,
        "axes.spines.right": False, "legend.frameon": False,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    })
# deterministic layout: no bbox_inches='tight' (anchored legends/footnote)
matplotlib.rcParams["savefig.bbox"] = "standard"
matplotlib.rcParams["pdf.fonttype"] = 42


def panel_tag(ax, tag, dx=-0.03, dy=1.03, fs=10.5):
    """Bold uppercase panel tag above the axes' top-left corner (matches the
    in-image (A)/(B) lettering used by main Figures 1-7)."""
    ax.text(dx, dy, tag, transform=ax.transAxes, fontsize=fs,
            fontweight="bold", ha="left", va="bottom")


# ---------- fixed naming / palettes (identical to Stage B deck)
COND_NOTE = "NEOLP n=3 (27,084 cells) vs EOLP n=2 (19,362 cells)"
SHORT = {0: "c0 Fibroblast SFRP2+", 1: "c1 Fibroblast CXCL14+", 2: "c2 Pericyte/SMC",
         3: "c3 T resting", 4: "c4 T/NK stressed", 5: "c5 T/NK low-complex.",
         6: "c6 CD8+ GZMK+ EM", 7: "c7 CD4+ Treg", 8: "c8 T naive IL7R+",
         9: "c9 Plasmablast GZMB+", 10: "c10 Macrophage LYZ+",
         11: "c11 Epith. basal/prolif.", 13: "c13 Mast", 14: "c14 B cell",
         15: "c15 Plasma cell", 16: "c16 CD8+ GZMB+/CXCL13+",
         17: "c17 Endo lymph.", 19: "c19 Endo blood", 20: "c20 Myeloid S100-hi"}
FAMILY_ORDER = ["T/NK", "B/Plasma", "Myeloid", "Mast", "Epithelial",
                "Fibroblast", "Pericyte/SMC", "Endothelial"]
FAM_COL = {"T/NK": "#882255", "B/Plasma": "#CC6677", "Myeloid": "#DDCC77",
           "Mast": "#117733", "Epithelial": "#88CCEE", "Fibroblast": "#44AA99",
           "Pericyte/SMC": "#999933", "Endothelial": "#332288"}
GENE_ORDER = ["GDF11", "SUOX", "EGFL8", "TAPBP", "HLA-DPA1"]
UP_C, DN_C = "#B02418", "#1F4E9C"      # EOLP up / down (panel a)
NS_C, ND_C = "#E9E9E9", "#8A8A8A"       # no supported change / not detected
SEG_C = {"comp": "#B0793A", "within": "#3E6FB0", "inter": "#9AA0A6"}

# ---------- data
g = pd.read_csv(os.path.join(TAB, "geo_gse211630_mr_backmap_gene.csv"))
st = pd.read_csv(os.path.join(TAB, "geo_gse211630_mr_backmap_subtype.csv"))
st["cluster"] = st["cluster"].astype(int)
st["supported"] = st["sig_cell"] | st["sig_perm"]
st["both"] = st["sig_cell"] & st["sig_perm"]
st["dir"] = np.sign(st["mean_all_EOLP"] - st["mean_all_NEOLP"])

# row order: family blocks, cluster ascending inside each block (as Fig B1)
cl2fam = dict(zip(st.cluster, st.family))
rows, ypos, seps = [], {}, []
yy = 0
for f in FAMILY_ORDER:
    cs = sorted(c for c in cl2fam if cl2fam[c] == f)
    for c in cs:
        ypos[int(c)] = yy
        rows.append(int(c))
        yy += 1
    seps.append(yy - 0.5)
nrows = yy - 1                       # 18

badge = {"GDF11": "risk \u2191 OR 2.88", "SUOX": "risk \u2191 OR 1.23",
         "EGFL8": "protective \u2193 OR 0.20", "TAPBP": "risk \u2191 OR 1.55",
         "HLA-DPA1": "protective \u2193 OR 0.88"}
badge_col = {"GDF11": "#505050", "SUOX": "#505050", "TAPBP": "#505050",
             "EGFL8": "#2E6B4F", "HLA-DPA1": "#2E6B4F"}

fig = plt.figure(figsize=(11.9, 7.5))
gs = fig.add_gridspec(1, 2, width_ratios=[1.18, 1.0],
                      left=0.125, right=0.985, top=0.855, bottom=0.26,
                      wspace=0.28)
axa = fig.add_subplot(gs[0])
axb = fig.add_subplot(gs[1])

# =====================================================================
# Panel (a): gene x subtype direction dot grid
# =====================================================================
xs = {gene: i for i, gene in enumerate(GENE_ORDER)}
for _, r in st.iterrows():
    x = xs[r["gene"]]; y = ypos[r["cluster"]]
    fr = r["expr_frac_pooled"]
    if r["n_expr_NEOLP"] + r["n_expr_EOLP"] == 0:
        axa.scatter(x, y, marker="x", s=26, color=ND_C, lw=0.9, zorder=3)
        continue
    if r["supported"]:
        face = UP_C if r["dir"] > 0 else (DN_C if r["dir"] < 0 else NS_C)
        ring = "black" if r["both"] else "white"
        rw = 0.9 if r["both"] else 0.35
    else:
        face, ring, rw = NS_C, "#B9B9B9", 0.35
    axa.scatter(x, y, s=3 + 170 * np.sqrt(fr), color=face,
                edgecolor=ring, linewidth=rw, zorder=3)

axa.set_xlim(-0.55, 4.55)
axa.set_ylim(-1.05, 20.6)
# gene headers above the top row
for i, gene in enumerate(GENE_ORDER):
    axa.text(i, 19.35, gene, ha="center", va="center",
             fontsize=7.6, fontweight="bold", style="italic")
    axa.text(i, 18.62, badge[gene], ha="center", va="center",
             fontsize=5.9, color=badge_col[gene])
axa.set_xticks([])
order_cl = [c for c, _p in sorted(ypos.items(), key=lambda kv: kv[1])]
axa.set_yticks(list(range(nrows + 1)))
axa.set_yticklabels([SHORT[c] for c in order_cl], fontsize=6.0)
# tint tick labels by family (as Fig B1)
for lab, c in zip(axa.get_yticklabels(), order_cl):
    fam = cl2fam[c]
    rgb = np.asarray(matplotlib.colors.to_rgb(FAM_COL[fam]))
    lab.set_color(matplotlib.colors.to_hex(np.clip(rgb * 0.55 + 0.05, 0, 1)))
axa.tick_params(axis="y", length=0, pad=1)
for sp in seps:
    axa.axhline(sp, color="#D6D6D6", lw=0.5, zorder=1)
axa.set_title("Target-message direction inside each cell subtype (EOLP vs NEOLP)",
              fontsize=8.2, loc="left", pad=2)
panel_tag(axa, "A", dx=-0.045, dy=1.05)
axa.set_ylabel("cell subtype (row colour = family, as in Fig B1)",
               fontsize=7.0, labelpad=6)

# legend A inside the lower-left whitespace of panel a
la = [
    Line2D([], [], marker="o", linestyle="none", ms=6.5, mfc=UP_C, mec="white",
           mew=0.5, label="EOLP up (supported)"),
    Line2D([], [], marker="o", linestyle="none", ms=6.5, mfc=DN_C, mec="white",
           mew=0.5, label="EOLP down (supported)"),
    Line2D([], [], marker="o", linestyle="none", ms=6.5, mfc=NS_C, mec="#B9B9B9",
           mew=0.5, label="no supported change"),
    Line2D([], [], marker="x", linestyle="none", ms=5.5, color=ND_C,
           label="not detected"),
    Line2D([], [], marker="o", linestyle="none", ms=6.5, mfc=NS_C, mec="black",
           mew=0.9, label="both evidence tracks"),
]
for fr in (0.05, 0.25, 0.75):
    la.append(Line2D([], [], marker="o", linestyle="none",
                     ms=2.4 * np.sqrt(3 + 170 * np.sqrt(fr)) * 0.62,
                     mfc="#B9B9B9", mec="none", label=f"{fr:.2f}"))
lgA = axa.legend(handles=la[:5], loc="upper left", bbox_to_anchor=(-0.02, -0.245),
                 fontsize=6.0, ncol=2, columnspacing=1.0, handletextpad=0.7,
                 labelspacing=0.45)
axa.add_artist(lgA)
axa.legend(handles=[la[5], la[6], la[7]], loc="upper left",
           bbox_to_anchor=(0.60, -0.245), fontsize=6.0, ncol=1,
           columnspacing=0.8, handletextpad=0.7, labelspacing=0.45,
           title="pooled expr. fraction", title_fontsize=6.0)

# =====================================================================
# Panel (b): per-gene horizontal waterfall (composition/within/inter)
# =====================================================================
yb = {gene: i for i, gene in enumerate(GENE_ORDER)}
terms = [("comp", "comp_contribution", "composition"),
         ("within", "within_contribution", "within-celltype"),
         ("inter", "interaction_contribution", "interaction")]
maxx = 0.0
for _, r in g.iterrows():
    cum = 0.0
    prev = 0.0
    for key, col, _lab in terms:
        v = r[col]
        cum += v
        bot, h = min(prev, cum), abs(v)
        axb.barh(yb[r["gene"]], h, left=bot,
                 color=SEG_C[key], height=0.55, edgecolor="none", zorder=2)
        prev = cum
    tot = r["tissue_delta"]
    axb.plot(tot, yb[r["gene"]], marker="D", ms=5.2, color="#111111",
             mec="white", mew=0.5, zorder=5)
    lbl = f"{r['tissue_rel_pct']:+.0f}%"
    # headline: relative tissue change (delta of diamond annotation)
    axb.text(tot, yb[r["gene"]] + 0.34, lbl, ha="center", va="bottom",
             fontsize=6.4, fontweight="bold", color="#111111")
    maxx = max(maxx, abs(cum))

axb.set_yticks(list(yb.values()))
axb.set_yticklabels([f"{gene}\n{badge[gene].split(' OR')[0]}" for gene in GENE_ORDER],
                    fontsize=6.4)
for lab in axb.get_yticklabels():
    lab.set_linespacing(1.0)
axb.set_ylim(-0.9, len(GENE_ORDER) - 0.35)
lim = 1.35 * maxx
axb.set_xlim(-lim * 0.6, lim)
axb.axvline(0, color="#999999", lw=0.6, zorder=1)
axb.set_xlabel(r"$\Delta$ tissue mean log1p expression  (EOLP $\minus$ NEOLP)",
               fontsize=7.0, labelpad=6)
axb.set_title("What moves the tissue-level target message: composition vs per-cell shift",
              fontsize=8.2, loc="left", pad=2)
panel_tag(axb, "B", dx=-0.045, dy=1.05)
axb.text(0.99, 1.03,
         "polarity assumption: circulating cis-pQTL direction extends to lesion "
         "mRNA (mRNA \u2260 protein; tissue \u2260 plasma)",
         transform=axb.transAxes, fontsize=5.7, color="#555555",
         ha="right", va="bottom")

lgB = [
    Patch(facecolor=SEG_C["comp"], label="composition (subtype fractions)"),
    Patch(facecolor=SEG_C["within"], label="within-celltype (per-cell expr.)"),
    Patch(facecolor=SEG_C["inter"], label="interaction"),
    Line2D([], [], marker="D", linestyle="none", ms=5, color="#111111",
           label="observed tissue \u0394"),
]
axb.legend(handles=lgB, loc="upper left", bbox_to_anchor=(-0.01, -0.26),
           fontsize=6.0, handletextpad=0.6, labelspacing=0.6)

fig.suptitle("Reverse-mapping of the five OLP MR drug-target genes onto the "
             "erosive-lesion cell atlas (GSE211630)",
             fontsize=9.4, weight="bold", x=0.01, ha="left", y=0.955)
fig.text(0.01, 0.920,
         "GSE211630 \u00b7 NEOLP n=3 (27,084 cells) vs EOLP n=2 (19,362 cells) "
         "\u00b7 tissue \u0394 = \u03a3 composition + \u03a3 within-celltype + "
         "\u03a3 interaction (Oaxaca, base = NEOLP) \u00b7 support = cell "
         "BH-FDR<0.05 OR sample-perm p\u22640.2; ring = both tracks.",
         fontsize=6.0, color="#555555")

out = os.path.join(FIGDIR, "geo_gse211630_mr_backmap")


# ---------- collision scan (canvas draw first; 2.5 px tolerance; own ticks OK)
def collision_scan(fig):
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
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


fig.savefig(out + ".png", dpi=300)
fig.savefig(out + ".pdf")
fig.savefig(os.path.join(QA, "fig_mr_backmap.png"), dpi=150)
collision_scan(fig)
plt.close(fig)
print("Fig C1 saved:", out + ".png/.pdf, QA 150 dpi")