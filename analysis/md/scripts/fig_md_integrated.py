#!/usr/bin/env python3
"""fig_md_integrated.py — the single docking+MD figure for the paper (nine panels, 175 mm).

Panels
  a  docking redone: published vs recomputed affinity, what dropped below the good cut-off
  b  reproducibility: seed-to-seed spread, and recovery of the crystallographic sulfate
  c  pose retention / protein flexibility / surface exposure (three stacked traces)
  d  the cofactor effect: same ligand, pocket emptied vs pocket with MTE+SO4
  e  contact hotspots: which residues each ligand actually touches
  f  departure, not rearrangement: ligand drift vs ligand-internal conformational RMSD
  g  hydrogen bonds: occupancy per complex over the trajectory
  h  MM-GBSA binding free energy per complex
  i  pocket druggability from fpocket, SUOX vs GDF11

data -> results/tables/*.csv + *.json, figures/md/figure_md_docking.{png,pdf}
"""
import csv
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

exec(compile(open(r"\\?\D:\wisp-science\skills\figure-style\runtime.py", encoding="utf-8").read(),
             "runtime.py", "exec"))
apply_figure_style(frame="boxed", sizes=(8.2, 7.0, 6.0))
BASE, LEG, TICK = 8.2, 7.0, 6.0          # the three-step ladder, named for readability below

# __file__ lives in <project>/analysis/md/scripts, so the project root is four levels up
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
T = lambda f: os.path.join(ROOT, "results", "tables", f)
OUT = os.path.join(ROOT, "figures", "main_figures_v2")
os.makedirs(OUT, exist_ok=True)

MM = 1 / 25.4
W, H = 175, 232

COMPLEXES = ["SUOX_Spiperone", "SUOX_Vorinostat", "SUOX_Sanguinarine", "GDF11_GW8510"]
COL = {"SUOX_Spiperone": "#1f77b4", "SUOX_Vorinostat": "#e8850c",
       "SUOX_Sanguinarine": "#17a2b8", "GDF11_GW8510": "#c0508f"}
LABEL = {"SUOX_Spiperone": "Spiperone", "SUOX_Vorinostat": "Vorinostat",
         "SUOX_Sanguinarine": "Sanguinarine", "GDF11_GW8510": "GW-8510"}

# ------------------------------------------------------------------ data
redo = list(csv.DictReader(open(T("docking_redo_comparison.csv"), encoding="utf-8")))
seed = json.load(open(os.path.join(ROOT, "analysis", "docking", "redo",
                                  "docking_redo_summary.json"), encoding="utf-8"))
so4 = json.load(open(os.path.join(ROOT, "analysis", "docking", "redo",
                                  "validation_so4.json"), encoding="utf-8"))
curves = np.load(T("md_curves.npy"), allow_pickle=True).item()
contacts = list(csv.DictReader(open(T("md_residue_contacts.csv"), encoding="utf-8")))
rmsf = list(csv.DictReader(open(T("md_rmsf.csv"), encoding="utf-8")))
sasa = list(csv.DictReader(open(T("md_sasa.csv"), encoding="utf-8")))
hbond = list(csv.DictReader(open(T("md_hbonds.csv"), encoding="utf-8")))
mmgbsa = list(csv.DictReader(open(T("md_mmgbsa.csv"), encoding="utf-8")))
fpx = json.load(open(T("fpocket_summary.json"), encoding="utf-8"))

def runs_of(cx):
    return [v for k, v in curves.items() if k.startswith(cx + "/")]

# ------------------------------------------------------------------ figure
fig = plt.figure(figsize=(W * MM, H * MM))
gs = fig.add_gridspec(3, 3, left=0.075, right=0.985, top=0.965, bottom=0.055,
                      wspace=0.42, hspace=0.52)

def tag(ax, letter, dx=-0.20, dy=1.14):
    ax.text(dx, dy, letter.upper(), transform=ax.transAxes, fontsize=13, fontweight="bold",
            ha="left", va="top")

# ============================================================ a  docking redone
ax = fig.add_subplot(gs[0, 0])
pub = np.array([float(r["published_kcal"]) for r in redo])
new = np.array([float(r["redo_best_kcal"]) for r in redo])
# "lost" = was under the -7 cut-off before the redo and is not any more
lost = np.array([int(r["published_pass"]) == 1 and int(r["redo_pass"]) == 0 for r in redo])
ax.scatter(pub[~lost], new[~lost], s=17, c="0.45", lw=0, label=f"unchanged ({(~lost).sum()})", zorder=3)
ax.scatter(pub[lost], new[lost], s=20, c="#d95f02", lw=0, label=f"lost good-docking ({lost.sum()})", zorder=4)
lim = (-9.3, -1.8)
ax.plot(lim, lim, ls="--", c="0.6", lw=0.8, zorder=1)
ax.axhline(-7, ls=":", c="0.55", lw=0.8, zorder=1)
ax.axvline(-7, ls=":", c="0.55", lw=0.8, zorder=1)
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("published affinity (kcal/mol)", fontsize=BASE, labelpad=1)
ax.set_ylabel("redone affinity (kcal/mol)", fontsize=BASE, labelpad=1)
ax.set_title(f"redo: {int((pub < -7).sum())} $\\rightarrow$ {int((new < -7).sum())} pass",
             fontsize=BASE, loc="left", pad=4)
ax.tick_params(labelsize=TICK)
ax.legend(fontsize=TICK, frameon=False, loc="upper left", handletextpad=0.4, borderaxespad=0.2)
tag(ax, "a")

# ============================================================ b  reproducibility
ax = fig.add_subplot(gs[0, 1])
labels, spreads = [], []
for tgt in ("SUOX", "GDF11"):
    spreads.append([r["seed_spread"] for r in seed if r["target"] == tgt])
    labels.append(tgt)
for gi, (tgt, vals) in enumerate(zip(labels, spreads)):
    y = np.full(len(vals), gi) + np.random.RandomState(0).uniform(-0.06, 0.06, len(vals))
    ax.scatter(vals, y, s=18, c=COL["SUOX_Spiperone" if gi == 0 else "GDF11_GW8510"], lw=0, zorder=3)
    ax.plot([np.median(vals)] * 2, [gi - 0.18, gi + 0.18], c="black", lw=1.4, zorder=4)
ax.set_yticks(range(len(labels)))
ax.set_yticklabels([f"{l} (n={len(v)})" for l, v in zip(labels, spreads)], fontsize=LEG)
ax.set_ylim(-0.6, len(labels) - 0.4)
ax.set_xlabel("seed-to-seed spread (kcal/mol)", fontsize=BASE, labelpad=1)
ax.set_title("redone docking is reproducible", fontsize=BASE, loc="left", pad=4)
ax.text(0.03, 0.97, f"SO$_4$ self-dock {so4['S_displacement_A']:.2f} Å (limit 2 Å)",
        transform=ax.transAxes, fontsize=TICK, color="0.35", ha="left", va="top")
ax.tick_params(labelsize=TICK)
tag(ax, "b")

# ============================================================ c  three stacked traces
gsc = gs[0, 2].subgridspec(3, 1, height_ratios=[1.0, 0.62, 0.62], hspace=0.30)

axc = fig.add_subplot(gsc[0])
for cx in COMPLEXES:
    for r in runs_of(cx):
        axc.plot(r["time"], r["drift"], color=COL[cx], lw=0.8, alpha=0.55, zorder=2)
axc.axhline(0.5, color="0.6", lw=0.8, ls="--", zorder=1)
axc.set_xlim(0, 200); axc.set_ylim(-0.2, 6.3)
axc.text(198, 0.58, "0.5 nm", fontsize=TICK, color="0.45", va="bottom", ha="right")
axc.set_ylabel("drift (nm)", fontsize=LEG, labelpad=1)
axc.set_title("only Spiperone holds its pose", fontsize=BASE, loc="left", pad=4)
axc.set_xticklabels([]); axc.tick_params(labelsize=TICK)
axc.legend(handles=[Line2D([], [], color=COL[c], lw=1.6, label=LABEL[c]) for c in COMPLEXES],
           fontsize=TICK - 0.6, frameon=False, loc="upper left", bbox_to_anchor=(1.015, 1.0),
           handletextpad=0.4, borderaxespad=0.0, labelspacing=0.32)

axr = fig.add_subplot(gsc[1])
for cx in COMPLEXES:
    per = {}
    for row in rmsf:
        if row["run"].startswith(cx + "/"):
            per.setdefault(int(row["resid"]), []).append(float(row["rmsf_nm"]))
    xs = sorted(per)
    axr.plot(xs, [np.mean(per[x]) for x in xs], color=COL[cx], lw=0.7, alpha=0.9)
axr.axhline(0.30, color="0.7", lw=0.7, ls=":", zorder=0)
axr.set_xlim(0, max(xs)); axr.set_ylim(0, 1.0)
axr.set_ylabel("RMSF (nm)", fontsize=LEG, labelpad=1)
axr.text(max(xs), 0.33, "0.3 nm", fontsize=TICK, color="0.55", ha="right", va="bottom")
axr.set_xticklabels([]); axr.tick_params(labelsize=TICK)

axs = fig.add_subplot(gsc[2])
for cx in COMPLEXES:
    per = {}
    for row in sasa:
        rn = row["run"]
        if rn.startswith(cx + "/") or rn.startswith(cx + "_"):
            per.setdefault(round(float(row["time_ns"])), []).append(float(row["sasa_protein_nm2"]))
    if not per:
        continue
    xs = sorted(per)
    axs.plot(xs, [np.mean(per[x]) for x in xs], color=COL[cx], lw=0.7, alpha=0.9)
axs.set_xlim(0, 200)
axs.set_xlabel("time (ns)", fontsize=BASE, labelpad=1)
axs.set_ylabel("SASA (nm$^2$)", fontsize=LEG, labelpad=1)
axs.tick_params(labelsize=TICK)
tag(axc, "c", dx=-0.30, dy=1.16)

# ============================================================ d  cofactor effect
ax = fig.add_subplot(gs[1, 0])
plain = {"SUOX_Sanguinarine": [4.233], "SUOX_Spiperone": [1.537], "SUOX_Vorinostat": [0.541]}
def end_drift(cx):
    return [r["drift"][-1] for r in runs_of(cx)]
tri = ["SUOX_Sanguinarine", "SUOX_Spiperone", "SUOX_Vorinostat"]
for i, cx in enumerate(tri):
    pf = plain[cx][0]
    ax.scatter([i], [pf], s=34, facecolors="none", edgecolors="#7aa6c2", lw=1.1, zorder=3)
    for v in end_drift(cx):
        ax.scatter([i + 0.22], [v], s=30, c=COL[cx], lw=0, zorder=4)
ax.axhline(0.5, ls="--", c="0.6", lw=0.8, zorder=1)
ax.set_xticks(range(len(tri)))
ax.set_xticklabels([LABEL[c] for c in tri], fontsize=LEG, rotation=18, ha="right",
                   rotation_mode="anchor")
ax.set_xlim(-0.5, len(tri) - 0.5)
ax.set_ylabel("drift at end of run (nm)", fontsize=BASE, labelpad=1)
ax.set_title("cofactor rescues Spiperone", fontsize=BASE, loc="left", pad=4)
ax.set_ylim(0, 5.6)
ax.legend(handles=[Line2D([], [], marker="o", ls="", mfc="none", mec="#7aa6c2",
                          label="no cofactor (60 ns)"),
                   Line2D([], [], marker="o", ls="", color="0.3",
                          label="MTE+SO$_4$ (200 ns)")],
          fontsize=TICK - 0.5, frameon=False, loc="upper left", handletextpad=0.2,
          borderaxespad=0.3, labelspacing=0.3)
ax.tick_params(labelsize=TICK)
tag(ax, "d")

# ============================================================ e  contact hotspots
ax = fig.add_subplot(gs[1, 1])
res = sorted({int(r["resid"]) for r in contacts})
keep = []
for rr in res:
    vals = [float(r["contact_freq"]) for r in contacts if int(r["resid"]) == rr]
    if max(vals) >= 0.30:
        keep.append(rr)
keep = keep[:10]
M = np.zeros((len(COMPLEXES), len(keep)))
for i, cx in enumerate(COMPLEXES):
    for j, rr in enumerate(keep):
        v = [float(r["contact_freq"]) for r in contacts
             if r["run"].startswith(cx + "/") and int(r["resid"]) == rr]
        M[i, j] = np.mean(v) if v else 0.0
im = ax.imshow(M, aspect="auto", cmap="viridis", vmin=0, vmax=1)
ax.set_yticks(range(len(COMPLEXES)))
ax.set_yticklabels([LABEL[c] for c in COMPLEXES], fontsize=LEG)
ax.set_xticks(range(len(keep)))
ax.set_xticklabels([str(r) for r in keep], fontsize=TICK, rotation=90)
ax.set_xlabel("SUOX residue", fontsize=BASE, labelpad=1)
ax.set_title("contact hotspots", fontsize=BASE, loc="left", pad=4)
cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
cb.ax.set_title("freq.", fontsize=TICK, pad=3)
cb.ax.tick_params(labelsize=TICK)
tag(ax, "e")

# ============================================================ f  departure vs rearrangement
ax = fig.add_subplot(gs[1, 2])
for cx in COMPLEXES:
    for r in runs_of(cx):
        ax.scatter(r["drift"], r["conf"], s=3, c=COL[cx], alpha=0.25, lw=0, zorder=3)
ax.axvline(0.5, ls="--", c="0.6", lw=0.8, zorder=1)
ax.set_xlabel("drift (nm)", fontsize=BASE, labelpad=1)
ax.set_ylabel("conformational RMSD (nm)", fontsize=BASE, labelpad=1)
ax.set_title("departure, not rearrangement", fontsize=BASE, loc="left", pad=4)
ax.text(0.62, 0.36, "pose lost", fontsize=TICK, color="0.4", ha="left")
ax.text(1.05, 0.06, "pose held", fontsize=TICK, color="0.4", ha="left")
ax.tick_params(labelsize=TICK)
tag(ax, "f")

# ============================================================ g  hydrogen bonds
ax = fig.add_subplot(gs[2, 0])
hb = np.zeros(len(COMPLEXES))
hb_sd = np.zeros(len(COMPLEXES))
for i, cx in enumerate(COMPLEXES):
    per_run = {}
    for row in hbond:
        if row["run"].startswith(cx + "/"):
            per_run.setdefault(row["run"], []).append(float(row["hbond_occ"]))
    vals = [np.sum(v) for v in per_run.values()]
    hb[i] = np.mean(vals) if vals else 0
    # panels G and H are labelled "mean +/- s.d."; np.std defaults to ddof=0
    # (population SD), which understates the sample SD by 0.82x at n = 3
    hb_sd[i] = np.std(vals, ddof=1) if len(vals) > 1 else 0
ax.bar(range(len(COMPLEXES)), hb, yerr=hb_sd, capsize=2.5,
       color=[COL[c] for c in COMPLEXES], edgecolor="white", lw=0.6, width=0.62)
ax.set_xticks(range(len(COMPLEXES)))
ax.set_xticklabels([LABEL[c] for c in COMPLEXES], fontsize=LEG, rotation=18, ha="right",
                   rotation_mode="anchor")
ax.set_ylabel("H-bonds per frame", fontsize=BASE, labelpad=1)
ax.set_title("polar contacts are few", fontsize=BASE, loc="left", pad=4)
ax.tick_params(labelsize=TICK)
tag(ax, "g")

# ============================================================ h  MM-GBSA
ax = fig.add_subplot(gs[2, 1])
mean, sdv = [], []
for cx in COMPLEXES:
    v = [float(r["dG_bind_kcal"]) for r in mmgbsa if r["complex"] == cx]
    mean.append(np.mean(v)); sdv.append(np.std(v, ddof=1) if len(v) > 1 else 0)
ypos = np.arange(len(COMPLEXES))[::-1]
ax.barh(ypos, mean, xerr=sdv, capsize=2.5, color=[COL[c] for c in COMPLEXES],
        edgecolor="white", lw=0.6, height=0.60)
ax.axvline(0, c="0.4", lw=0.8)
ax.set_yticks(ypos)
ax.set_yticklabels([LABEL[c] for c in COMPLEXES], fontsize=LEG)
ax.set_xlabel("MM-GBSA $\\Delta G_{bind}$ (kcal/mol)", fontsize=BASE, labelpad=1)
ax.set_title("binding energy does not track stability", fontsize=BASE - 0.6, loc="left", pad=4)
ax.tick_params(labelsize=TICK)
# the number goes outside the error bar; only the longest bar carries it inside, in white,
# so it cannot collide with the y tick labels
for y, m, sd in zip(ypos, mean, sdv):
    if m < -12:
        ax.text(m * 0.45, y, f"{m:.1f}", fontsize=TICK, ha="center", va="center", color="white")
    else:
        ax.text(m - sd - 0.45, y, f"{m:.1f}", fontsize=TICK, ha="right", va="center", color="0.15")
ax.set_xlim(-19.5, 0.5)
tag(ax, "h")

# ============================================================ i  pocket druggability
ax = fig.add_subplot(gs[2, 2])
names = ["SUOX", "GDF11"]
drug = [fpx[n]["top"]["druggability"] for n in names]
vol = [fpx[n]["top"]["volume"] for n in names]
bars = ax.bar(range(2), drug, color=[COL["SUOX_Spiperone"], COL["GDF11_GW8510"]],
              edgecolor="white", lw=0.6, width=0.55)
for i, (b, d, v, n) in enumerate(zip(bars, drug, vol, names)):
    ax.text(i, d + 0.03, f"{d:.2f}", ha="center", fontsize=TICK, color="0.2")
    ax.text(i, 0.04, f"{v:.0f} Å$^3$\n{fpx[n]['n_pockets']} pockets", ha="center",
            fontsize=TICK - 0.4, color="white")
ax.set_xticks(range(2)); ax.set_xticklabels(names, fontsize=LEG)
ax.set_ylim(0, 1.12)
ax.set_ylabel("druggability score", fontsize=BASE, labelpad=1)
ax.set_title("SUOX has the ligandable pocket", fontsize=BASE, loc="left", pad=4)
ax.tick_params(labelsize=TICK)
tag(ax, "i")

png = os.path.join(OUT, "fig04_docking_md.png")
fig.savefig(png, dpi=300, bbox_inches="tight", pad_inches=0.02)
fig.savefig(png.replace(".png", ".pdf"), bbox_inches="tight", pad_inches=0.02)
print(f"wrote {os.path.relpath(png, ROOT)}")

# ---- overlap check (crude but catches the usual offenders)
fig.canvas.draw()
r = fig.canvas.get_renderer()
boxes = []
for a in fig.axes:
    for t in [a.title] + list(a.texts) + [a.xaxis.label, a.yaxis.label]:
        if t is None or not t.get_text().strip():
            continue
        try:
            bb = t.get_window_extent(r)
        except Exception:
            continue
        if bb.width > 0 and bb.height > 0:
            boxes.append((t.get_text().strip()[:22], bb))
hits = []
for i in range(len(boxes)):
    for j in range(i + 1, len(boxes)):
        if boxes[i][1].overlaps(boxes[j][1]):
            hits.append((boxes[i][0], boxes[j][0]))
print(f"text/text overlaps: {len(hits)}")
for h in hits[:10]:
    print("   ", h)
