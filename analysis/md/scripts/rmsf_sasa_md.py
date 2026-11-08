#!/usr/bin/env python3
"""rmsf_sasa_md.py — per-residue RMSF over the production runs.

Completes panel c of the integrated figure, which currently shows only ligand drift.

RMSF is taken from the Cα atoms after superposing every frame on frame 0, at 1 ns spacing (201
frames per 200 ns run). SASA is NOT computed here: this MDAnalysis build ships without the SASA
module, so it is produced by `gmx sasa` on the server (60_sasa_gmx.sh) and merged by merge_sasa.py.

Outputs
  results/tables/md_rmsf.csv       run, resid, resname, rmsf_nm
"""
import csv
import glob
import os
import sys

import numpy as np
import MDAnalysis as mda

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
BASE = os.path.join(ROOT, "remote", "easyai", "md_solute_20261018", "unpacked")
STRIDE = 10          # 1 ns at a 100 ps trajectory


def runs():
    out = []
    for top in sorted(glob.glob(os.path.join(BASE, "complexes", "*", "*", "solute_top.gro"))):
        d = os.path.dirname(top)
        xtc = os.path.join(d, "solute.xtc")
        if os.path.exists(xtc) and os.path.getsize(xtc) > 100_000:
            out.append((os.path.basename(os.path.dirname(d)), os.path.basename(d), top, xtc))
    return out


def analyse(top, xtc):
    u = mda.Universe(top, xtc)
    prot = u.select_atoms("protein")
    lig = u.select_atoms("resname LIG or resname UNL")
    if len(lig) == 0:
        return None
    ca = u.select_atoms("protein and name CA")

    # --- RMSF: superpose each frame on frame 0, then track each Cα ---
    u.trajectory[0]
    ref = ca.positions.copy()
    devs = []
    sasa_complex, sasa_prot, sasa_lig, times = [], [], [], []
    for ts in u.trajectory[::STRIDE]:
        # Kabsch superposition of the Cα set onto frame 0, then the per-residue deviation.
        # (MDAnalysis's own aligner wants two Universes; a direct Kabsch is simpler and exact.)
        P = ca.positions
        Pc = P - P.mean(0); Rc = ref - ref.mean(0)
        V, S, Wt = np.linalg.svd(Pc.T @ Rc)
        d = np.sign(np.linalg.det(V @ Wt))
        Rot = V @ np.diag([1.0, 1.0, d]) @ Wt
        Pa = Pc @ Rot
        devs.append(np.linalg.norm(Pa - Rc, axis=1))
        times.append(ts.time / 1000.0)
        # SASA is not computed here: this MDAnalysis build has no SASA module, so it comes from
        # `gmx sasa` (see 60_sasa_gmx.sh) and is merged by merge_sasa.py.
        pass

    devs = np.asarray(devs)                # (nframes, nres)
    rmsf = np.sqrt((devs ** 2).mean(0)) / 10.0   # A -> nm
    return dict(ca_resids=[int(r) for r in ca.residues.resids],
                ca_names=list(ca.residues.resnames),
                rmsf=rmsf,
                times=np.asarray(times),
                sasa_complex=np.asarray(sasa_complex),
                sasa_prot=np.asarray(sasa_prot))


def main():
    rs = runs()
    print(f"runs: {len(rs)}")
    rmsf_rows = []
    for cx, run, top, xtc in rs:
        r = analyse(top, xtc)
        if r is None:
            print(f"  {cx}/{run}: no ligand, skipped")
            continue
        label = f"{cx}/{run}"
        for rid, nm, val in zip(r["ca_resids"], r["ca_names"], r["rmsf"]):
            rmsf_rows.append(dict(run=label, resid=rid, resname=nm, rmsf_nm=round(float(val), 4)))
        # SASA is reported in A^2 by MDAnalysis; convert to nm^2 and report the buried surface
        # SASA rows come from merge_sasa.py (gmx sasa), not from here
        print(f"  {label:44s} RMSF mean {r['rmsf'].mean():.3f} nm  "
              f"max {r['rmsf'].max():.3f} nm")

    os.makedirs(os.path.join(ROOT, "results", "tables"), exist_ok=True)
    p1 = os.path.join(ROOT, "results", "tables", "md_rmsf.csv")
    with open(p1, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["run", "resid", "resname", "rmsf_nm"])
        w.writeheader(); w.writerows(rmsf_rows)
    print(f"\nwrote {os.path.relpath(p1, ROOT)} ({len(rmsf_rows)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
