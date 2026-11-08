# -*- coding: utf-8 -*-
"""
Fig D1 — GSE211630 OLP single-cell expansion, Stage D (LIANA/CellChat) main figure.
Cells: cell-cell communication (liana rank_aggregate, consensus resource) between
NEOLP vs EOLP, restricted to curated focus receivers (T/NK clusters 3-8,16 plus
the two fibroblast hosts c0/c1 that carry Stage-C risk-gene loss).

Canonical classes (from geo_gse211630_liana_signal_shift.csv, column 'class'):
  EOLP_lost   = significant only in NEOLP
  sustained   = significant in BOTH conditions (also both-present by construction)
  EOLP_gained = significant only in EOLP
  none        = not significant in either

Panels
  a  presence-flux heatmap: incoming focus L-R pairs per receiver x class (+ net flux)
  b  strength shift: sum(delta_prob) over sustained pairs per receiver (EOLP - NEOLP)
  c  ligand extremes: ligand-level sum(delta_prob) over sustained pairs targeting
     the 7 T/NK focus receivers

Outputs
  figures/geo_gse211630_liana_comm.png|.pdf          (final)
  analysis/scrna/_qa/figD1_{full,a,b,c}.png          (QA crops)
  results/tables/geo_gse211630_liana_figd1_summary.csv  (canonical numbers for QA/README)

Dependencies (pixi env analysis/scrna): python, pandas, numpy, matplotlib.
figure-style sidecar is loaded when reachable for house style; otherwise a
minimal equivalent style is applied so the script remains standalone.
"""
from __future__ import annotations
import sys, os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[2]
TAB = ROOT / "results" / "tables"
FIG = ROOT / "figures"
QA = ROOT / "analysis" / "scrna" / "_qa"
QA.mkdir(parents=True, exist_ok=True)

# ----------------------------------------------------------------- house style
_SKILL_RT = r"\\?\D:\wisp-science\skills\figure-style\runtime.py"
_has_style = False
if os.path.exists(_SKILL_RT):
    try:
        exec(compile(open(_SKILL_RT, encoding="utf-8").read(), "fs_runtime.py", "exec"), globals())
        apply_figure_style(frame="open", sizes=(8, 7, 6))
        _has_style = True
    except Exception as exc:  # pragma: no cover - fall back below
        print(f"[figD1] figure-style sidecar load failed ({exc}); using fallback style")
if not _has_style:  # minimal equivalent fallback
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
        "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.6, "xtick.direction": "out", "ytick.direction": "out",
        "axes.titleweight": "regular",
    })

# ----------------------------------------------------------------- canonical inputs
sh = pd.read_csv(TAB / "geo_gse211630_liana_signal_shift.csv", encoding="utf-8-sig")
inf = pd.read_csv(TAB / "geo_gse211630_liana_incoming_focus.csv", encoding="utf-8-sig")
ann = pd.read_csv(TAB / "geo_gse211630_celltype_annotation.csv", encoding="utf-8-sig")

TN_FOCUS = {3, 4, 5, 6, 7, 8, 16}          # T/NK focus receivers (Stage B subtypes)
RECEIVERS = [16, 6, 1, 3, 8, 4, 7, 5, 0]   # row order: descending sum(delta_prob)

_SUB = {r: (ann.loc[ann["leiden_0.5"] == r, "celltype_sub"].iloc[0]) for r in RECEIVERS}
_FAM = {r: (ann.loc[ann["leiden_0.5"] == r, "celltype_family"].iloc[0]) for r in RECEIVERS}
SHORT = {
    0:  "c0  Fib SFRP2+",
    1:  "c1  Fib CXCL14+",
    3:  "c3  T resting CD4/CD8",
    4:  "c4  T/NK HSP-stressed",
    5:  "c5  T/NK low-complexity",
    6:  "c6  CD8 GZMK+ EM",
    7:  "c7  CD4 Treg",
    8:  "c8  T naive-like IL7Rhi",
    16: "c16 CD8 GZMB+/CXCL13+ exh",
}

# ---- a) presence flux (n = L-R pairs, summed over source families) ----------
pres = (inf.groupby(["target_cluster", "class"])["n"].sum().unstack(fill_value=0)
        .reindex(columns=["EOLP_lost", "sustained", "EOLP_gained", "none"], fill_value=0)
        .loc[RECEIVERS])
pres["net_flux"] = pres["EOLP_gained"] - pres["EOLP_lost"]

# ---- b) strength shift over sustained pairs ---------------------------------
foc = sh[sh["receiver_focus"] & (sh["class"] == "sustained")].copy()
stren = (foc.groupby("target_cluster")
         .agg(n_pairs=("delta_prob", "size"), sum_dp=("delta_prob", "sum"))
         .loc[RECEIVERS])
stren["is_tnk"] = stren.index.isin(TN_FOCUS)

# ---- c) ligand extremes, sustained pairs targeting T/NK focus receivers -----
lt = (foc[foc["target_cluster"].isin(TN_FOCUS)]
      .groupby("ligand")
      .agg(n_pairs=("delta_prob", "size"), sum_dp=("delta_prob", "sum"))
      .query("n_pairs >= 14")
      .sort_values("sum_dp"))
lig_bottom = lt.head(6).iloc[::-1]          # strongest declines
lig_top = lt.tail(4)                        # strongest gains
lig_show = pd.concat([lig_bottom, lig_top])
# consolidate the three classical MHC-I ligands shown
mhc1_sum = lt.loc[[g for g in ["HLA-A", "HLA-C", "HLA-E"] if g in lt.index], "sum_dp"].sum()

# canonical summary for QA / README (long tail of extremes beyond figure)
summ = {"receiver": RECEIVERS}
figd1 = pd.DataFrame({
    "receiver": RECEIVERS,
    "family": [_FAM[r] for r in RECEIVERS],
    "subtype": [_SUB[r] for r in RECEIVERS],
    "focus_TNK": [r in TN_FOCUS for r in RECEIVERS],
    "pairs_EOLP_lost": pres["EOLP_lost"].values,
    "pairs_sustained": pres["sustained"].values,
    "pairs_EOLP_gained": pres["EOLP_gained"].values,
    "pairs_none": pres["none"].values,
    "net_flux": pres["net_flux"].values,
    "n_sustained_sig_pairs": stren["n_pairs"].values,
    "sum_delta_prob_sustained": stren["sum_dp"].values,
})
lig_df = pd.DataFrame({
    "ligand": lig_show.index,
    "n_sustained_pairs_TNK": lig_show["n_pairs"].values,
    "sum_delta_prob_sustained_TNK": lig_show["sum_dp"].values,
})
lig_df.to_csv(TAB / "geo_gse211630_liana_figd1_summary_ligands.csv",
              index=False, encoding="utf-8-sig")
figd1.to_csv(TAB / "geo_gse211630_liana_figd1_summary.csv",
             index=False, encoding="utf-8-sig")
print("[figD1] summary CSVs written")
print(figd1[["receiver", "focus_TNK", "net_flux", "n_sustained_sig_pairs",
             "sum_delta_prob_sustained"]].to_string(index=False, float_format=lambda v: f"{v:.3f}"))
print("MHC-I (HLA-A/C/E) sustained sum_dp toward T/NK focus =", f"{mhc1_sum:.3f}")

# ----------------------------------------------------------------- colors / style constants
NEG = "#377eb8"     # decline (EOLP < NEOLP)
POS = "#ff7f00"     # gain
GRAY = "#b8b8b8"
CMAP = matplotlib.colormaps["YlGnBu"]

# ----------------------------------------------------------------- figure
W_IN, H_IN = 7.15, 5.9
fig = plt.figure(figsize=(W_IN, H_IN), dpi=300)
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.14], width_ratios=[1.0, 0.92],
                      left=0.155, right=0.88, top=0.90, bottom=0.15,
                      hspace=0.66, wspace=0.5)

# ---------------- panel a: presence-flux heatmap ------------------------------
axA = fig.add_subplot(gs[0, :])
cls_cols = ["EOLP_lost", "sustained", "EOLP_gained"]
cls_lbl = {"EOLP_lost": "lost in EOLP", "sustained": "sustained",
           "EOLP_gained": "gained in EOLP"}
X = pres[cls_cols].values.astype(float)
vmax = X.max()
im = axA.imshow(X, cmap=CMAP, norm=Normalize(0, vmax), aspect="auto")
axA.set_xticks(range(4))
axA.set_xticklabels([cls_lbl[c] for c in cls_cols] + ["net flux"])
axA.set_yticks(range(len(RECEIVERS)))
axA.set_yticklabels([SHORT[r] for r in RECEIVERS])
for i in range(X.shape[0]):
    for j in range(3):
        v = X[i, j]
        col = "white" if v > 0.55 * vmax else "black"
        axA.text(j, i, f"{int(v)}", ha="center", va="center", fontsize=6, color=col)
# net-flux column: light cells + sign-coloured values (magnitude printed)
for i, r in enumerate(RECEIVERS):
    nf = int(pres.loc[r, "net_flux"])
    axA.add_patch(matplotlib.patches.Rectangle(
        (2.5, i - 0.5), 1.0, 1.0, facecolor="#f2f2f2", edgecolor="#d8d8d8", lw=0.5))
    col = POS if nf > 0 else (NEG if nf < 0 else "0.45")
    axA.text(3.0, i, f"{nf:+d}", ha="center", va="center", fontsize=6.5,
             color=col, fontweight="bold")
axA.axvline(2.5, color="0.85", lw=0.8)
axA.set_xlim(-0.5, 3.5)
axA.grid(False)
axA.set_title("Incoming focus signalling per receiver by presence class (n = L–R pairs)",
              loc="left")
for sp in axA.spines.values():
    sp.set_visible(True)
axA.spines["top"].set_visible(False); axA.spines["right"].set_visible(False)
cbar = fig.colorbar(im, ax=axA, fraction=0.016, pad=0.02)
cbar.set_label("n pairs (counts)", fontsize=7)
cbar.ax.tick_params(labelsize=6)
cbar.set_ticks(np.linspace(0, vmax, 5))
cbar.set_ticklabels([f"{int(t)}" for t in np.linspace(0, vmax, 5)])
panel_letter(axA, "a", case="upper") if "panel_letter" in globals() else None

# ---------------- panel b: strength shift (sustained pairs) -------------------
axB = fig.add_subplot(gs[1, 0])
y = np.arange(len(RECEIVERS))[::-1]           # same order top->bottom as panel a
vals = stren["sum_dp"].values
for yi, r, v in zip(y, RECEIVERS, vals):
    is_tnk = r in TN_FOCUS
    col = NEG if v < 0 else POS
    if is_tnk:
        face, edge = col, col
    else:
        face = "#e2ebf5" if v < 0 else "#fdecdd"
        edge = GRAY
    axB.barh(yi, v, height=0.62, color=face, edgecolor=edge, linewidth=0.8,
             hatch="" if is_tnk else "//")
    axB.text(v + np.sign(v) * 0.18, yi, f"{v:+.2f}", va="center",
             ha="left" if v > 0 else "right", fontsize=6, color=col)
axB.axvline(0, color="k", lw=0.7)
axB.set_yticks(y)
axB.set_yticklabels([SHORT[r] for r in RECEIVERS], fontsize=6)
axB.set_xlim(-8.7, 1.5)
axB.set_xlabel("Σ Δinteraction probability, sustained pairs (EOLP − NEOLP)")
axB.set_title("Sustained T/NK-directed signalling declines in EOLP",
              loc="left")
panel_letter(axB, "b", case="upper") if "panel_letter" in globals() else None

# ---------------- panel c: ligand extremes ------------------------------------
axC = fig.add_subplot(gs[1, 1])
yv = np.arange(len(lig_show))
for yi, (lig, row) in enumerate(lig_show.iterrows()):
    v = row["sum_dp"]
    col = NEG if v < 0 else POS
    axC.barh(yi, v, height=0.6, color=col)
    axC.text(v + np.sign(v) * 0.09, yi, f"{v:+.2f}", va="center",
             ha="left" if v > 0 else "right", fontsize=6, color=col)
axC.axvline(0, color="k", lw=0.7)
axC.set_yticks(yv)
axC.set_yticklabels([f"{g}   (n={int(lig_show.loc[g, 'n_pairs'])})" for g in lig_show.index],
                    fontsize=6)
axC.set_xlim(-3.0, 0.7)
axC.set_xlabel("Σ Δprob per ligand (sustained → T/NK focus receivers)")
axC.set_title("Sustained losses concentrated in MHC-I ligands",
              loc="left")
panel_letter(axC, "c", case="upper") if "panel_letter" in globals() else None

# ---------------- shared sign legend (bottom) ---------------------------------
leg = [Patch(facecolor=NEG, label="loss / decline in EOLP"),
       Patch(facecolor=POS, label="gain in EOLP"),
       Patch(facecolor="#f2f2f2", edgecolor=GRAY, hatch="//",
             label="hatched bars: non-T/NK receivers c0/c1 (B)")]
fig.legend(handles=leg, loc="lower center", bbox_to_anchor=(0.5, 0.012),
           frameon=False, ncol=3, fontsize=6, handlelength=1.6,
           handleheight=0.9, columnspacing=1.8, labelspacing=0.4)

# ----------------------------------------------------------------- save + QA
out_png = FIG / "geo_gse211630_liana_comm.png"
out_pdf = FIG / "geo_gse211630_liana_comm.pdf"
fig.savefig(out_png, dpi=300)
fig.savefig(out_pdf)
print("[figD1] saved", out_png, out_pdf)

# 1) collision scan (draw first so window extents share the final layout space)
fig.canvas.draw()
rend = fig.canvas.get_renderer()
labels = [(t, t.get_window_extent(rend)) for t in fig.findobj(matplotlib.text.Text)
          if t.get_text().strip() and t.get_visible()]
frames = [(s, s.get_window_extent(rend)) for axx in fig.axes
          for s in axx.spines.values() if s.get_visible()]
own_ticks = {axx: set(axx.get_xticklabels(which="both") + axx.get_yticklabels(which="both"))
             for axx in fig.axes}
hits = [(a, b) for i, (a, ba) in enumerate(labels)
        for b, bb in labels[i + 1:] if ba.overlaps(bb)]
hits += [(t, s) for t, bt in labels for s, bs in frames
         if bt.overlaps(bs) and t not in own_ticks[s.axes]]
if hits:
    def _txt(x):
        return x.get_text() if hasattr(x, "get_text") else "<spine>"
    for a, b in hits[:20]:
        ba = a.get_window_extent(rend) if hasattr(a, "get_window_extent") else None
        bb = b.get_window_extent(rend) if hasattr(b, "get_window_extent") else None
        ea = (f"bbox=({ba.x0:.0f},{ba.y0:.0f},{ba.x1:.0f},{ba.y1:.0f})" if ba else "")
        eb = (f"bbox=({bb.x0:.0f},{bb.y0:.0f},{bb.x1:.0f},{bb.y1:.0f})" if bb else "")
        print("[figD1] COLLISION:", repr(_txt(a))[:46], ea, "<->", repr(_txt(b))[:46], eb)
    raise SystemExit(f"[figD1] collision scan failed: {len(hits)} hits")
print("[figD1] collision scan clean")

# 2) panel crops for visual QA (pixel crops of the saved PNG)
from PIL import Image as _PIL  # noqa: E402
_img = _PIL.open(out_png).convert("RGB")
for axx, nm in [(axA, "figD1_a.png"), (axB, "figD1_b.png"), (axC, "figD1_c.png")]:
    bb = axx.get_tightbbox(rend).expanded(1.14, 1.1)
    H = _img.height
    x0, y0 = max(0, int(bb.x0) - 6), max(0, H - int(bb.y1) - 6)
    x1, y1 = min(_img.width, int(bb.x1) + 6), min(H, H - int(bb.y0) + 6)
    _img.crop((x0, y0, x1, y1)).save(QA / nm)
    print("[figD1] crop", nm)
import shutil as _sh
_sh.copyfile(out_png, QA / "figD1_full.png")
print("[figD1] QA crops written to", QA)
