# -*- coding: utf-8 -*-
"""Back-map + LIANA plate, split into two files.

Main Figure 7 (was "A/B" of the merged plate): original Fig11 (MR back-map).
  A target-message direction dot grid | B Oaxaca tissue decomposition
  -> figures/main_figures/fig07_backmap_liana.{png,pdf}  (12.2 x 6.6 in)

Supplementary Figure 10 (was C/D/E of the merged plate): original Fig12 (LIANA).
  A presence-flux heatmap (old C) | B sustained-pair strength shift (old D)
  | C ligand extremes (old E)
  -> figures/supplementary/supp_fig10_liana_communication.{png,pdf} (12.2 x 9.4 in)

Renumbered 2026-09-13: this plate moved from Supplementary Figure 7 to
Supplementary Figure 10 when the supplementary set was reordered to follow the
manuscript narrative; the plate content is unchanged.

Restructured 2026-09-30 for the main-figure re-plan: each panel is now a function that draws into
a caller-supplied Axes, so panels can be recomposed into new main figures while both existing
plates stay unchanged. `main()` still creates the original gridspecs and passes the original
rectangles/letters, so it reproduces all four output files exactly.

Both plates share one data-load block. Reads the same Stage-C / Stage-D tables and
reuses the exact drawing code from plot_mr_backmap_gse211630.py and
scrna_figD1_build_gse211630.py (verbatim, re-tagged).
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
os.chdir(ROOT)
TAB = os.path.join("results", "tables")
FIGDIR = os.path.join("figures", "main_figures")      # main plate (A, B)
SUPPDIR = os.path.join("figures", "supplementary")    # supplementary plate (old C-E)
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(SUPPDIR, exist_ok=True)
QA = os.path.join("analysis", "scrna", "_qa")
os.makedirs(QA, exist_ok=True)

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7,
                     "axes.linewidth": 0.6, "svg.fonttype": "none"})

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
UP_C, DN_C = "#B02418", "#1F4E9C"      # EOLP up / down (panel A)
NS_C, ND_C = "#E9E9E9", "#8A8A8A"       # no supported change / not detected
SEG_C = {"comp": "#B0793A", "within": "#3E6FB0", "inter": "#9AA0A6"}
NEG = "#377eb8"                          # LIANA decline (EOLP < NEOLP)
POS = "#ff7f00"                          # LIANA gain
GRAY = "#b8b8b8"
CMAP = matplotlib.colormaps["YlGnBu"]
TN_FOCUS = {3, 4, 5, 6, 7, 8, 16}
RECEIVERS = [16, 6, 1, 3, 8, 4, 7, 5, 0]

# WIDE odds ratios from the corrected allele alignment. Higher genetically predicted
# expression goes with LOWER OLP risk for GDF11, SUOX and TAPBP, and with HIGHER risk
# for EGFL8 and HLA-DPA1.
badge = {"GDF11": "protective \u2193 OR 0.35", "SUOX": "protective \u2193 OR 0.81",
         "EGFL8": "risk \u2191 OR 4.95", "TAPBP": "protective \u2193 OR 0.65",
         "HLA-DPA1": "risk \u2191 OR 1.14"}
badge_col = {"GDF11": "#2E6B4F", "SUOX": "#2E6B4F", "TAPBP": "#2E6B4F",
             "EGFL8": "#505050", "HLA-DPA1": "#505050"}

SHORT_D = {
    0:  "c0  Fib SFRP2+",
    1:  "c1  Fib CXCL14+",
    3:  "c3  T resting CD4/CD8",
    4:  "c4  T/NK HSP-stressed",
    5:  "c5  T/NK low-complex.",
    6:  "c6  CD8 GZMK+ EM",
    7:  "c7  CD4 Treg",
    8:  "c8  T naive IL7Rhi",
    16: "c16 CD8 GZMB+/CXCL13+",
}

# Original plate geometry (inches). main() reuses these so both published plates are reproduced
# exactly; the recomposition builder supplies its own figures and rectangles.
PLATE_MAIN_IN = (12.2, 6.6)
GS_MAIN = dict(width_ratios=[1.02, 1.0], wspace=0.30,
               left=0.14, right=0.90, top=0.95, bottom=0.26)
PLATE_SUPP_IN = (12.2, 9.4)
GS_SUPP = dict(height_ratios=[2.50, 3.05], hspace=0.60,
               left=0.14, right=0.90, top=0.93, bottom=0.155)
GS_SUPP_ROW2 = dict(width_ratios=[1.0, 1.0], wspace=0.50)

SUPP_NOTE = ("A\u2013C: LIANA rank_aggregate consensus; sustained = significant in both conditions; "
             "T/NK focus receivers = c3\u2013c8 + c16.")


def panel_tag(ax, tag, dx=-0.02, dy=1.02, fs=11):
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


# ------------------------------------------------------------------ data
g = pd.read_csv(os.path.join(TAB, "geo_gse211630_mr_backmap_gene.csv"))
st = pd.read_csv(os.path.join(TAB, "geo_gse211630_mr_backmap_subtype.csv"))
st["cluster"] = st["cluster"].astype(int)
st["supported"] = st["sig_cell"] | st["sig_perm"]
st["both"] = st["sig_cell"] & st["sig_perm"]
st["dir"] = np.sign(st["mean_all_EOLP"] - st["mean_all_NEOLP"])

sh = pd.read_csv(os.path.join(TAB, "geo_gse211630_liana_signal_shift.csv"), encoding="utf-8-sig")
inf = pd.read_csv(os.path.join(TAB, "geo_gse211630_liana_incoming_focus.csv"), encoding="utf-8-sig")
ann = pd.read_csv(os.path.join(TAB, "geo_gse211630_celltype_annotation.csv"), encoding="utf-8-sig")

# panel C data
pres = (inf.groupby(["target_cluster", "class"])["n"].sum().unstack(fill_value=0)
        .reindex(columns=["EOLP_lost", "sustained", "EOLP_gained", "none"], fill_value=0)
        .loc[RECEIVERS])
pres["net_flux"] = pres["EOLP_gained"] - pres["EOLP_lost"]

# panel D data
foc = sh[sh["receiver_focus"] & (sh["class"] == "sustained")].copy()
stren = (foc.groupby("target_cluster")
         .agg(n_pairs=("delta_prob", "size"), sum_dp=("delta_prob", "sum"))
         .loc[RECEIVERS])
stren["is_tnk"] = stren.index.isin(TN_FOCUS)

# panel E data
lt = (foc[foc["target_cluster"].isin(TN_FOCUS)]
      .groupby("ligand")
      .agg(n_pairs=("delta_prob", "size"), sum_dp=("delta_prob", "sum"))
      .query("n_pairs >= 14")
      .sort_values("sum_dp"))
lig_bottom = lt.head(6).iloc[::-1]
lig_top = lt.tail(4)
lig_show = pd.concat([lig_bottom, lig_top])

# back-map row order (family blocks)
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
nrows = yy - 1


# =====================================================================
# Panel A: gene x subtype direction dot grid
# =====================================================================
def panel_A(ax, tag="A", badge_fontsize=5.9, gene_fontsize=7.6, ylabel_fontsize=6.0,
            size_key=True, legend_bbox=(-0.02, -0.16),
            legend_fontsize=6.0, legend_labelspacing=0.45, legend_ncol=2,
            title="Target-message direction inside each cell subtype (EOLP vs NEOLP)"):
    """Target-message direction dot grid.

    badge_fontsize / gene_fontsize / ylabel_fontsize let a recomposition shrink the column headers
    and row labels when the panel is narrower than its original half of a 12.2 in plate, and
    size_key=False drops the dot-area key (which sits to the right, below the panel, where a
    recomposition may need that space). The defaults reproduce the original figure.
    """
    xs = {gene: i for i, gene in enumerate(GENE_ORDER)}
    for _, r in st.iterrows():
        x = xs[r["gene"]]; y = ypos[r["cluster"]]
        fr = r["expr_frac_pooled"]
        if r["n_expr_NEOLP"] + r["n_expr_EOLP"] == 0:
            ax.scatter(x, y, marker="x", s=26, color=ND_C, lw=0.9, zorder=3)
            continue
        if r["supported"]:
            face = UP_C if r["dir"] > 0 else (DN_C if r["dir"] < 0 else NS_C)
            ring = "black" if r["both"] else "white"
            rw = 0.9 if r["both"] else 0.35
        else:
            face, ring, rw = NS_C, "#B9B9B9", 0.35
        ax.scatter(x, y, s=3 + 170 * np.sqrt(fr), color=face,
                   edgecolor=ring, linewidth=rw, zorder=3)

    ax.set_xlim(-0.55, 4.55)
    ax.set_ylim(-1.05, 20.6)
    for i, gene in enumerate(GENE_ORDER):
        ax.text(i, 19.35, gene, ha="center", va="center",
                fontsize=gene_fontsize, fontweight="bold", style="italic")
        ax.text(i, 18.62, badge[gene], ha="center", va="center",
                fontsize=badge_fontsize, color=badge_col[gene])
    ax.set_xticks([])
    order_cl = [c for c, _p in sorted(ypos.items(), key=lambda kv: kv[1])]
    ax.set_yticks(list(range(nrows + 1)))
    ax.set_yticklabels([SHORT[c] for c in order_cl], fontsize=ylabel_fontsize)
    for lab, c in zip(ax.get_yticklabels(), order_cl):
        fam = cl2fam[c]
        rgb = np.asarray(matplotlib.colors.to_rgb(FAM_COL[fam]))
        lab.set_color(matplotlib.colors.to_hex(np.clip(rgb * 0.55 + 0.05, 0, 1)))
    ax.tick_params(axis="y", length=0, pad=1)
    for sp in seps:
        ax.axhline(sp, color="#D6D6D6", lw=0.5, zorder=1)
    ax.set_title(title, fontsize=8.2, loc="left", pad=2)
    if tag:
        panel_tag(ax, tag, dx=-0.05, dy=1.03)
    ax.set_ylabel("cell subtype (row colour = family)",
                  fontsize=7.0, labelpad=6)

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
    # Size key: marker areas mirror panel A (s = 3 + 170*sqrt(fr)); the legend draws
    # diameters, so ms = sqrt(s) * 1.00 -- the key circles are exactly the plotted
    # sizes. The former 1.49x factor made the largest circle (18.2 pt) exceed the
    # legend row pitch, so the three key circles merged; labelspacing 1.25 now keeps
    # clear white gaps (25 px and 13 px at 300 dpi).
    for fr in (0.05, 0.25, 0.75):
        la.append(Line2D([], [], marker="o", linestyle="none",
                         ms=np.sqrt(3 + 170 * np.sqrt(fr)) * 1.00,
                         mfc="#B9B9B9", mec="none", label=f"{fr:.2f}"))
    # legend_bbox is a parameter because the two-row key hangs below the axes; the band gap under
    # this panel decides how close to the axes it has to sit.
    lgA = ax.legend(handles=la[:5], loc="upper left", bbox_to_anchor=legend_bbox,
                    ncol=legend_ncol, columnspacing=1.0 if legend_ncol > 1 else 2.0,
                    fontsize=legend_fontsize, handletextpad=0.7,
                    labelspacing=legend_labelspacing)
    ax.add_artist(lgA)
    if size_key:
        ax.legend(handles=[la[5], la[6], la[7]], loc="upper left",
                  bbox_to_anchor=(0.62, -0.16), fontsize=6.0, ncol=1,
                  columnspacing=0.8, handletextpad=0.7, labelspacing=1.25, borderpad=0.6,
                  title="pooled expr. fraction", title_fontsize=6.0)
    return ax


# =====================================================================
# Panel B: per-gene horizontal waterfall
# =====================================================================
def panel_B(ax, tag="B", title="Tissue target message: composition vs per-cell shift",
            legend_loc="upper left", legend_bbox=(-0.01, -0.20), legend_fontsize=6.0):
    """Oaxaca tissue decomposition waterfall; `title` can be shortened for a narrower panel."""
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
            ax.barh(yb[r["gene"]], h, left=bot,
                    color=SEG_C[key], height=0.55, edgecolor="none", zorder=2)
            prev = cum
        tot = r["tissue_delta"]
        ax.plot(tot, yb[r["gene"]], marker="D", ms=5.2, color="#111111",
                mec="white", mew=0.5, zorder=5)
        lbl = f"{r['tissue_rel_pct']:+.0f}%"
        ax.text(tot, yb[r["gene"]] + 0.34, lbl, ha="center", va="bottom",
                fontsize=6.4, fontweight="bold", color="#111111")
        maxx = max(maxx, abs(cum))

    ax.set_yticks(list(yb.values()))
    ax.set_yticklabels([f"{gene}\n{badge[gene].split(' OR')[0]}" for gene in GENE_ORDER],
                       fontsize=6.4)
    for lab in ax.get_yticklabels():
        lab.set_linespacing(1.0)
    ax.set_ylim(-0.9, len(GENE_ORDER) - 0.35)
    lim = 1.35 * maxx
    ax.set_xlim(-lim * 0.6, lim)
    ax.axvline(0, color="#999999", lw=0.6, zorder=1)
    ax.set_xlabel(r"$\Delta$ tissue mean log1p expression  (EOLP $\minus$ NEOLP)",
                  fontsize=7.0, labelpad=6)
    ax.set_title(title,
                 fontsize=8.2, loc="left", pad=2)
    if tag:
        panel_tag(ax, tag, dx=-0.05, dy=1.03)

    lgB = [
        Patch(facecolor=SEG_C["comp"], label="composition (subtype fractions)"),
        Patch(facecolor=SEG_C["within"], label="within-celltype (per-cell expr.)"),
        Patch(facecolor=SEG_C["inter"], label="interaction"),
        Line2D([], [], marker="D", linestyle="none", ms=5, color="#111111",
               label="observed tissue \u0394"),
    ]
    # legend_loc / legend_fontsize are parameters because a recomposition may need the key inside
    # the axes: parked below the axes it ran into the band underneath. Rows 2-5 extend left
    # (TAPBP -11%, EGFL8 -14%, SUOX -35%, GDF11 -18%), so the lower-right of this panel is empty.
    ax.legend(handles=lgB, loc=legend_loc, bbox_to_anchor=legend_bbox,
              fontsize=legend_fontsize, handletextpad=0.6, labelspacing=0.6)
    return ax


# =====================================================================
# Panel C: presence-flux heatmap (full width)
# =====================================================================
def panel_C(ax, tag="A", colorbar=True,
            title="Incoming focus signalling per receiver by presence class (n = L\u2013R pairs)",
            cbar_ax=None):
    """Presence-flux heatmap; `title` can be shortened for a narrower panel.

    Pass `cbar_ax` to place the colour bar in its own axes. Without it, matplotlib reserves the
    bar's space *inside* `ax`, which shortens this panel relative to its row neighbours.
    """
    cls_cols = ["EOLP_lost", "sustained", "EOLP_gained"]
    cls_lbl = {"EOLP_lost": "lost in EOLP", "sustained": "sustained",
               "EOLP_gained": "gained in EOLP"}
    X = pres[cls_cols].values.astype(float)
    vmax = X.max()
    im = ax.imshow(X, cmap=CMAP, norm=Normalize(0, vmax), aspect="auto")
    ax.set_xticks(range(4))
    ax.set_xticklabels([cls_lbl[c] for c in cls_cols] + ["net flux"])
    ax.set_yticks(range(len(RECEIVERS)))
    ax.set_yticklabels([SHORT_D[r] for r in RECEIVERS])
    for i in range(X.shape[0]):
        for j in range(3):
            v = X[i, j]
            col = "white" if v > 0.55 * vmax else "black"
            ax.text(j, i, f"{int(v)}", ha="center", va="center", fontsize=6, color=col)
    for i, r in enumerate(RECEIVERS):
        nf = int(pres.loc[r, "net_flux"])
        ax.add_patch(matplotlib.patches.Rectangle(
            (2.5, i - 0.5), 1.0, 1.0, facecolor="#f2f2f2", edgecolor="#d8d8d8", lw=0.5))
        col = POS if nf > 0 else (NEG if nf < 0 else "0.45")
        ax.text(3.0, i, f"{nf:+d}", ha="center", va="center", fontsize=6.5,
                color=col, fontweight="bold")
    ax.axvline(2.5, color="0.85", lw=0.8)
    ax.set_xlim(-0.5, 3.5)
    ax.grid(False)
    ax.set_title(title, loc="left", fontsize=8.2, pad=2)
    for sp in ax.spines.values():
        sp.set_visible(True)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

    cbar = None
    if colorbar:
        if cbar_ax is not None:
            cbar = ax.figure.colorbar(im, cax=cbar_ax)
        else:
            cbar = ax.figure.colorbar(im, ax=ax, fraction=0.016, pad=0.02)
        cbar.set_label("n pairs (counts)", fontsize=7)
        cbar.ax.tick_params(labelsize=6)
        cbar.set_ticks(np.linspace(0, vmax, 5))
        cbar.set_ticklabels([f"{int(t)}" for t in np.linspace(0, vmax, 5)])
    if tag:
        panel_tag(ax, tag, dx=-0.05, dy=1.03)
    return dict(ax=ax, cbar=cbar, im=im)


# =====================================================================
# Panel D: strength shift (sustained pairs)
# =====================================================================
def panel_D(ax, tag="B", title="Sustained T/NK-directed signalling declines in EOLP"):
    """Sustained-pair strength shift; `title` can be shortened for a narrower panel."""
    yD = np.arange(len(RECEIVERS))[::-1]
    vals = stren["sum_dp"].values
    for yi, r, v in zip(yD, RECEIVERS, vals):
        is_tnk = r in TN_FOCUS
        col = NEG if v < 0 else POS
        if is_tnk:
            face, edge = col, col
        else:
            face = "#e2ebf5" if v < 0 else "#fdecdd"
            edge = GRAY
        ax.barh(yi, v, height=0.62, color=face, edgecolor=edge, linewidth=0.8,
                hatch="" if is_tnk else "//")
        ax.text(v + np.sign(v) * 0.18, yi, f"{v:+.2f}", va="center",
                ha="left" if v > 0 else "right", fontsize=6, color=col)
    ax.axvline(0, color="k", lw=0.7)
    ax.set_yticks(yD)
    ax.set_yticklabels([SHORT_D[r] for r in RECEIVERS], fontsize=6)
    ax.set_xlim(-8.7, 1.5)
    ax.set_xlabel("\u03a3 \u0394interaction probability, sustained pairs (EOLP \u2212 NEOLP)")
    ax.set_title(title, loc="left", fontsize=8.2, pad=2)
    if tag:
        panel_tag(ax, tag, dx=-0.05, dy=1.03)
    return ax


# =====================================================================
# Panel E: ligand extremes
# =====================================================================
def panel_E(ax, tag="C"):
    yv = np.arange(len(lig_show))
    for yi, (lig, row) in enumerate(lig_show.iterrows()):
        v = row["sum_dp"]
        col = NEG if v < 0 else POS
        ax.barh(yi, v, height=0.6, color=col)
        ax.text(v + np.sign(v) * 0.09, yi, f"{v:+.2f}", va="center",
                ha="left" if v > 0 else "right", fontsize=6, color=col)
    ax.axvline(0, color="k", lw=0.7)
    ax.set_yticks(yv)
    ax.set_yticklabels([f"{g}   (n={int(lig_show.loc[g, 'n_pairs'])})" for g in lig_show.index],
                       fontsize=6)
    ax.set_xlim(-3.0, 0.7)
    ax.set_xlabel("\u03a3 \u0394prob per ligand (sustained \u2192 T/NK focus receivers)")
    ax.set_title("Sustained losses concentrated in MHC-I ligands",
                 loc="left", fontsize=8.2, pad=2)
    if tag:
        panel_tag(ax, tag, dx=-0.05, dy=1.03)
    return ax


def supp_plate_legend(figS):
    """The D/E colour key belongs to the supplementary plate (it explains its B and C)."""
    leg = [Patch(facecolor=NEG, label="loss / decline in EOLP"),
           Patch(facecolor=POS, label="gain in EOLP"),
           Patch(facecolor="#f2f2f2", edgecolor=GRAY, hatch="//",
                 label="hatched bars: non-T/NK receivers c0/c1 (B)")]
    return figS.legend(handles=leg, loc="lower center", bbox_to_anchor=(0.5, 0.008),
                       frameon=False, ncol=3, fontsize=6, handlelength=1.6,
                       handleheight=0.9, columnspacing=1.8, labelspacing=0.4)


def main():
    # Main plate: panels A and B only (C-E moved to the supplementary plate).
    # Panel geometry (axes widths, wspace, per-row labels) is the same as before;
    # only the height of the plate and the band below the axes changed, so the
    # two legends anchored under A and B keep their previous offsets.
    fig = plt.figure(figsize=PLATE_MAIN_IN)
    gs1 = fig.add_gridspec(1, 2, **GS_MAIN)
    panel_A(fig.add_subplot(gs1[0]), tag="A")
    panel_B(fig.add_subplot(gs1[1]), tag="B")

    # Supplementary plate: old C -> A (full width), old D -> B, old E -> C.
    figS = plt.figure(figsize=PLATE_SUPP_IN)
    gsS = figS.add_gridspec(2, 1, **GS_SUPP)
    panel_C(figS.add_subplot(gsS[0]), tag="A")
    gs3 = gsS[1].subgridspec(1, 2, **GS_SUPP_ROW2)
    panel_D(figS.add_subplot(gs3[0]), tag="B")
    panel_E(figS.add_subplot(gs3[1]), tag="C")

    supp_plate_legend(figS)

    # (no figure-level title on the main plate: the caption carries the heading and
    #  the polarity assumption; the grey note line was removed on request)
    figS.text(0.02, 0.975, SUPP_NOTE, fontsize=6.0, color="#555555")

    fig.savefig(os.path.join(FIGDIR, "fig07_backmap_liana.png"), dpi=300)
    fig.savefig(os.path.join(FIGDIR, "fig07_backmap_liana.pdf"))
    fig.savefig(os.path.join(QA, "fig7_main_merged.png"), dpi=150)
    figS.savefig(os.path.join(SUPPDIR, "supp_fig10_liana_communication.png"), dpi=300)
    figS.savefig(os.path.join(SUPPDIR, "supp_fig10_liana_communication.pdf"))
    figS.savefig(os.path.join(QA, "supp_fig10_liana.png"), dpi=150)
    collision_scan(fig, fig.canvas.get_renderer())
    collision_scan(figS, figS.canvas.get_renderer())
    plt.close("all")
    print("Main Fig7 (A/B)  -> figures/main_figures/fig07_backmap_liana.png")
    print("Supp Fig10 (A-C) -> figures/supplementary/supp_fig10_liana_communication.png")
    return fig, figS


if __name__ == "__main__":
    main()
