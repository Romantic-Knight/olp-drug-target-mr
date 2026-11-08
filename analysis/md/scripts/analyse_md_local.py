#!/usr/bin/env python3
"""analyse_md_local.py — stability curves, PCA and FEL from the harvested trajectories.

Everything here runs on the local copy under remote/easyai/md_solute_20261018/, so it does not need
the rental server. The trajectories are already PBC-corrected and contain only protein+ligand
(solute_top.gro gives the matching topology).

Outputs
  results/tables/md_stability_all.csv        one row per run, summary metrics
  results/tables/md_pca_summary.csv          PCA explained variance + ligand-region sampling
  figures/md/stability_curves.png            drift / conf / contacts vs time, coloured by condition
  figures/md/pca_scatter.png                 PC1-PC2 projections per complex
  figures/md/fel_<complex>.png               free-energy landscape from PC1/PC2
"""
import glob
import os
import sys

import numpy as np
import MDAnalysis as mda

# __file__ is analysis/md/scripts/<this file>, so the project root is three levels up
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
BASE = os.path.join(ROOT, "remote", "easyai", "md_solute_20261018", "unpacked")
TS = os.path.join(ROOT, "remote", "easyai", "md_solute_20261018", "analysis", "md_timeseries.csv")


def kabsch_rmsd(P, Q):
    """RMSD after optimal superposition of P onto Q (both Nx3)."""
    Pc = P - P.mean(0)
    Qc = Q - Q.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    R = V @ D @ Wt
    Prot = Pc @ R
    return float(np.sqrt(((Prot - Qc) ** 2).sum(1).mean()))


def collect_runs():
    """Every run directory that has both a solute trajectory and its matching topology."""
    runs = []
    for top in sorted(glob.glob(os.path.join(BASE, "complexes", "*", "*", "solute_top.gro"))):
        d = os.path.dirname(top)
        xtc = os.path.join(d, "solute.xtc")
        if os.path.exists(xtc) and os.path.getsize(xtc) > 100000:
            runs.append((os.path.basename(os.path.dirname(d)),
                         os.path.basename(d), top, xtc))
    return runs


def analyse_run(cx, run, top, xtc, stride=10):
    u = mda.Universe(top, xtc)
    # Two different selections on purpose: the backbone is what we superpose on, but contacts must
    # be counted against ALL protein atoms. Using the CA selection for both made every run look
    # non-contacting (4 atoms instead of 59) — CA atoms are rarely within 4 A of a buried ligand.
    bb = u.select_atoms("protein and name CA")
    prot = u.select_atoms("protein")
    lig = u.select_atoms("resname LIG or resname UNL")
    if len(lig) == 0:
        # the ligand may have been written under its own residue name; take the last residue
        lig = u.select_atoms(f"resid {u.residues[-1].resid}")
    if len(bb) == 0 or len(lig) == 0:
        return None

    ref_bb = bb.positions.copy()
    ref_lig = lig.positions.copy()

    times, drift, conf, cont = [], [], [], []
    for ts in u.trajectory[::stride]:
        b = bb.positions
        l = lig.positions
        # MDAnalysis works in Angstrom; the whole project reports nm
        drift.append(_rmsd_after_fit(b, l, ref_bb, ref_lig) / 10.0)
        conf.append(kabsch_rmsd(l, ref_lig) / 10.0)
        d = np.linalg.norm(prot.positions[:, None, :] - l[None, :, :], axis=2)
        cont.append(int((d < 4.0).any(axis=1).sum()))
        times.append(ts.time / 1000.0)

    out = dict(complex=cx, run=run, ns=float(times[-1]),
               time=np.array(times), drift=np.array(drift),
               conf=np.array(conf), contacts=np.array(cont))
    return out


def _rmsd_after_fit(P, L, Pref, Lref):
    """RMSD of the ligand after superposing the protein backbone onto its reference."""
    Pc = P - P.mean(0)
    Prefc = Pref - Pref.mean(0)
    V, S, Wt = np.linalg.svd(Pc.T @ Prefc)
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1.0, 1.0, d]) @ Wt
    Lm = (L - P.mean(0)) @ R
    Lrefc = Lref - Pref.mean(0)
    return float(np.sqrt(((Lm - Lrefc) ** 2).sum(1).mean()))


def pca_run(top, xtc, stride=10):
    """PCA on backbone CA coordinates; returns PC1/PC2 projections and explained variance."""
    u = mda.Universe(top, xtc)
    ca = u.select_atoms("protein and name CA")
    n = len(ca)
    X = []
    for _ in u.trajectory[::stride]:
        X.append(ca.positions.flatten())
    X = np.asarray(X)
    X -= X.mean(0)
    # SVD on the covariance (n_frames x 3N) is cheaper this way round
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    var = (S ** 2) / (S ** 2).sum()
    proj = U[:, :2] * S[:2]
    return proj, var[:5]


def fel(proj, bins=40):
    """Histogram-based free energy from PC1/PC2, in kT."""
    H, xe, ye = np.histogram2d(proj[:, 0], proj[:, 1], bins=bins)
    P = H / H.sum()
    with np.errstate(divide="ignore"):
        G = -np.log(P)
    G -= np.nanmin(G[np.isfinite(G)])
    return G, xe, ye


def main():
    runs = collect_runs()
    print(f"runs with a matching topology: {len(runs)}")
    if not runs:
        print("nothing to analyse - check the unpacked directory")
        return 1

    os.makedirs(os.path.join(ROOT, "results", "tables"), exist_ok=True)
    os.makedirs(os.path.join(ROOT, "figures", "md"), exist_ok=True)

    rows, curves, pcas = [], {}, {}
    for cx, run, top, xtc in runs:
        r = analyse_run(cx, run, top, xtc)
        if r is None:
            print(f"  {cx}/{run}: no ligand found, skipped")
            continue
        key = f"{cx}/{run}"
        curves[key] = r
        rows.append(dict(complex=cx, run=run, ns=round(r["ns"], 1),
                         drift_last=round(float(r["drift"][-1]), 3),
                         drift_mean=round(float(r["drift"].mean()), 3),
                         conf_last=round(float(r["conf"][-1]), 3),
                         contacts_last=int(r["contacts"][-1]),
                         contacts_mean=int(r["contacts"].mean())))
        try:
            proj, var = pca_run(top, xtc)
            pcas[cx] = (proj, var)
        except Exception as e:
            print(f"  {key}: PCA skipped ({type(e).__name__})")
        print(f"  {key:44s} {r['ns']:6.1f} ns  drift {r['drift'][-1]:.3f}  contacts {r['contacts'][-1]}")

    import csv
    with open(os.path.join(ROOT, "results", "tables", "md_stability_all.csv"), "w",
              newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote results/tables/md_stability_all.csv ({len(rows)} runs)")

    np.save(os.path.join(ROOT, "results", "tables", "md_curves.npy"),
            {k: v for k, v in curves.items()}, allow_pickle=True)
    print("wrote results/tables/md_curves.npy (raw curves for plotting)")

    # ---- plots ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # 1. stability curves
    # One colour per run, reused across the three panels and fed to the legend explicitly.
    # Plotting without capturing the Line2D handles let matplotlib's own colour cycle drift out of
    # sync with the labels, so the legend named colours the curves did not have.
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6))
    palette = plt.get_cmap("tab20").colors
    handles, labels = [], []
    for i, (key, r) in enumerate(sorted(curves.items())):
        col = palette[i % len(palette)]
        short = key.replace("SUOX_", "").replace("GDF11_", "GDF11 ")
        h, = axes[0].plot(r["time"], r["drift"], lw=1.0, color=col, label=short)
        axes[1].plot(r["time"], r["conf"], lw=1.0, color=col)
        axes[2].plot(r["time"], r["contacts"], lw=1.0, color=col)
        handles.append(h); labels.append(short)
    for ax, ttl, yl in zip(axes,
                           ("ligand drift (fit on protein)", "ligand conformational RMSD",
                            "protein-ligand contacts < 0.4 nm"),
                           ("RMSD (nm)", "RMSD (nm)", "count")):
        ax.set_title(ttl, fontsize=8)
        ax.set_xlabel("time (ns)", fontsize=8)
        ax.set_ylabel(yl, fontsize=8)
        ax.tick_params(labelsize=7)
        ax.axhline(0.5 if "drift" in ttl else 0, color="0.8", lw=0.8, ls="--", zorder=0)
    axes[0].legend(handles, labels, fontsize=4.6, frameon=True, framealpha=0.85,
                   ncol=2, handlelength=1.2, columnspacing=0.8, loc="upper left")
    axes[0].set_ylim(bottom=-0.15)
    fig.tight_layout()
    fig.savefig(os.path.join(ROOT, "figures", "md", "stability_curves.png"), dpi=300)
    plt.close(fig)
    print("wrote figures/md/stability_curves.png")

    # 2. PCA scatter + FEL per complex
    for cx, (proj, var) in sorted(pcas.items()):
        fig, ax = plt.subplots(figsize=(3.2, 2.8))
        ax.scatter(proj[:, 0], proj[:, 1], s=3, alpha=0.5, lw=0, c=np.arange(len(proj)), cmap="viridis")
        ax.set_xlabel(f"PC1 ({var[0]*100:.0f}%)", fontsize=8)
        ax.set_ylabel(f"PC2 ({var[1]*100:.0f}%)", fontsize=8)
        ax.set_title(cx, fontsize=8)
        ax.tick_params(labelsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(ROOT, "figures", "md", f"pca_{cx}.png"), dpi=300)
        plt.close(fig)

        G, xe, ye = fel(proj)
        fig, ax = plt.subplots(figsize=(3.2, 2.8))
        im = ax.pcolormesh(xe, ye, G.T, cmap="viridis", shading="auto", vmin=0, vmax=6)
        cb = fig.colorbar(im, ax=ax, pad=0.02)
        cb.set_label("free energy (kT)", fontsize=7)
        cb.ax.tick_params(labelsize=6)
        ax.set_xlabel(f"PC1 ({var[0]*100:.0f}%)", fontsize=8)
        ax.set_ylabel(f"PC2 ({var[1]*100:.0f}%)", fontsize=8)
        ax.set_title(f"FEL — {cx}", fontsize=8)
        ax.tick_params(labelsize=7)
        fig.tight_layout()
        fig.savefig(os.path.join(ROOT, "figures", "md", f"fel_{cx}.png"), dpi=300)
        plt.close(fig)
        print(f"wrote figures/md/pca_{cx}.png and fel_{cx}.png")

    return 0


if __name__ == "__main__":
    sys.exit(main())
