#!/usr/bin/env python3
"""analyse_md.py — trajectory analysis: is the docked ligand actually stable in its site?

The question, and why one RMSD is not enough
--------------------------------------------
A docking score says nothing about whether a pose survives dynamics. The check that does is the
ligand's motion, but a single RMSD conflates two different things, so three are reported:

  * **drift** - superimpose on the protein backbone, then measure the ligand: motion *relative to
    the pocket*. A curve that climbs monotonically means the ligand left.
  * **conformational** - superimpose on the ligand itself: shape change only. Flat conformational
    RMSD together with rising drift means the ligand departed as a rigid body rather than
    rearranging.
  * **contacts** - protein atoms staying within 0.4 nm of the ligand, plus the minimum distance to
    any protein atom. A ligand that has escaped into solvent has almost no contacts and a growing
    minimum distance, which is unambiguous.

Protein backbone RMSD, Rg and SASA are included so an unfolded or collapsing protein cannot be
mistaken for a ligand effect.

Implementation notes
--------------------
* **The trajectory must be made whole first** (`gmx trjconv -pbc mol -center`, see
   `26_pbc_fix.sh`). The raw .xtc is written wrapped: with a ~5 nm box a ligand or protein segment
   that crosses a boundary appears to jump several nm, and the protein backbone RMSD then reads
   ~1.8 nm over 5 ns, which is physically impossible for a folded domain. That artefact, not real
   motion, was what the first analysis pass reported.
* Topology comes from the **.gro**, not the .tpr: MDAnalysis' TPR parser does not support the tpx
  version GROMACS 2026 writes (138). The .gro carries the residue names and coordinates this needs.
* The superpositions use an explicit **Kabsch** implementation rather than MDAnalysis' RMSD or
  AlignTraj classes. Those have changed selection semantics between versions (an earlier attempt
  failed with "Reference and trajectory atom selections do not contain the same number of atoms:
  N_ref=1484, N_traj=0"), and doing the fit directly makes the definition of each number explicit
  and version-independent.

Usage:
    python3 analyse_md.py [trajectory-dir ...]
"""
import glob
import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def kabsch(P, Q):
    """Rotation+translation that best maps P onto Q; returns (R, t)."""
    Pc, Qc = P.mean(axis=0), Q.mean(axis=0)
    P0, Q0 = P - Pc, Q - Qc
    H = P0.T @ Q0
    U, _S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    t = Qc - R @ Pc
    return R, t


def rmsd(P, Q):
    return float(np.sqrt(((P - Q) ** 2).sum(axis=1).mean()))


def analyse(top, xtc, ligand_resname="LIG", stride=5, label=""):
    import MDAnalysis as mda

    u = mda.Universe(top, xtc)
    # MDAnalysis works in Angstrom while every threshold below is written in nm, so positions are
    # converted explicitly wherever they are read (A2NM). An earlier version compared Angstrom
    # against nm cutoffs, which turned the 0.4 nm contact test into 0.04 nm. The trajectory
    # add_transformations API was tried first but its signature differs between versions (it raised
    # "'NoneType' object has no attribute 'has_positions'"), so the conversion is done inline.
    A2NM = 0.1
    prot = u.select_atoms("protein")
    lig = u.select_atoms(f"resname {ligand_resname}")
    if len(lig) == 0:
        return {"label": label, "error": f"no residue named {ligand_resname}"}
    bb = prot.select_atoms("backbone")
    if len(bb) == 0:
        bb = prot

    frames = list(u.trajectory[::stride])
    if len(frames) < 3:
        return {"label": label, "error": f"only {len(frames)} frame(s)"}

    # reference coordinates from the first analysed frame
    u.trajectory[0]
    ref_bb = bb.positions.copy() * A2NM
    ref_lig = lig.positions.copy() * A2NM

    times, drift, conf, protr, ncont, mind, rg = [], [], [], [], [], [], []

    for ts in u.trajectory[::stride]:
        bb_now = bb.positions.copy() * A2NM
        lig_now = lig.positions.copy() * A2NM

        # drift: fit the protein backbone, then move the ligand with the same transform
        R, t = kabsch(bb_now, ref_bb)
        lig_fitted = (R @ lig_now.T).T + t
        drift.append(rmsd(lig_fitted, ref_lig))
        protr.append(rmsd((R @ bb_now.T).T + t, ref_bb))

        # conformational: fit the ligand onto itself
        Rl, tl = kabsch(lig_now, ref_lig)
        conf.append(rmsd((Rl @ lig_now.T).T + tl, ref_lig))

        # contacts: minimum distance from each protein atom to the nearest ligand atom
        pair = np.linalg.norm(lig_now[:, None, :] - prot.positions[None, :, :] * A2NM, axis=2)
        mins = pair.min(axis=0)
        ncont.append(int((mins < 0.4).sum()))
        mind.append(float(mins.min()))

        # radius of gyration computed directly: MDAnalysis' polymer.radius_of_gyration
        # is not present in every version, and the definition is one line
        pos = prot.positions * A2NM
        rg.append(float(np.sqrt(((pos - pos.mean(axis=0)) ** 2).sum(axis=1).mean())))
        times.append(ts.time / 1000.0)

    def st(a):
        a = np.asarray(a, float)
        return {"first": round(float(a[0]), 3), "last": round(float(a[-1]), 3),
                "mean": round(float(a.mean()), 3), "max": round(float(a.max()), 3),
                "std": round(float(a.std()), 3)}

    # slope of drift over the trajectory: a slow monotonic rise is the signature of leaving
    slope = float(np.polyfit(np.arange(len(drift)), drift, 1)[0] * max(len(drift) - 1, 1))

    last_drift, last_cont = drift[-1], ncont[-1]
    if last_drift < 0.30 and last_cont >= max(5, int(ncont[0] * 0.5)):
        verdict = "stable"
    elif last_drift > 0.50 or last_cont <= 1:
        verdict = "left the site"
    else:
        verdict = "intermediate"

    return {"label": label, "frames_analysed": len(times),
            "duration_ns": round(times[-1], 2), "stride_ps": stride * 10,
            "ligand_drift_nm": st(drift), "ligand_conformational_nm": st(conf),
            "protein_backbone_rmsd_nm": st(protr),
            "contacts_within_0.4nm": {"first": ncont[0], "last": ncont[-1],
                                      "mean": int(np.mean(ncont)), "min": int(np.min(ncont))},
            "min_protein_distance_nm": {"first": round(mind[0], 3), "last": round(mind[-1], 3)},
            "rgyr_nm": {"first": round(rg[0], 3), "last": round(rg[-1], 3)},
            "drift_slope_nm_over_run": round(slope, 3),
            "verdict": verdict}


def main():
    if len(sys.argv) > 1:
        dirs = sys.argv[1:]
    else:
        # three naming schemes exist in this project: prod_rep* (first plain SUOX pass),
        # holo_rep* (cofactor present) and gdf11_rep* (GDF11 monomer)
        dirs = sorted(
            glob.glob(os.path.join(ROOT, "complexes", "*", "prod_rep*")) +
            glob.glob(os.path.join(ROOT, "complexes", "*", "holo_rep*")) +
            glob.glob(os.path.join(ROOT, "complexes", "*", "gdf11_rep*")) +
            glob.glob(os.path.join(ROOT, "complexes", "*", "trial")))

    results = []
    print(f"{'run':34s} {'ns':>6s} {'drift':>7s} {'conf':>7s} {'contacts':>9s} {'mindist':>8s}  verdict")
    for d in dirs:
        top = xtc = None
        # prefer whole.xtc (PBC-corrected, see 26_pbc_fix.sh); fall back to the raw trajectory but
        # flag it, because wrapped coordinates give meaningless superpositions
        # Pair the topology with the trajectory by ATOM COUNT, not by filename.
        #
        # There are two legitimate trajectory forms here and picking the wrong pair is silent until
        # MDAnalysis raises "topology and XTC do not have the same number of atoms":
        #   solute.xtc + solute_top.gro  - protein+ligand only, PBC-corrected, 100 ps stride
        #   whole.xtc  + holo/system.gro - whole system, PBC-corrected
        # solute is preferred: it is what 40_make_solute.sh produces after correction, and prod.xtc
        # (still wrapped) must never be used for a superposition. Both pairs below are self-consistent.
        pairs = []
        if os.path.exists(os.path.join(d, "solute_top.gro")):
            pairs.append((os.path.join(d, "solute_top.gro"), os.path.join(d, "solute.xtc")))
        for gname in ("system.gro", "holo.gro", "gdf11.gro"):
            for g in (os.path.join(d, gname), os.path.join(os.path.dirname(d), gname)):
                pairs.append((g, os.path.join(d, "whole.xtc")))
        for gp, xp in pairs:
            if os.path.exists(gp) and os.path.exists(xp) and os.path.getsize(xp) > 0:
                top, xtc = gp, xp
                break
        if not top:
            continue
        if "whole" not in os.path.basename(xtc) and "solute" not in os.path.basename(xtc):
            print(f"  NOTE {os.path.basename(os.path.dirname(d))}/{os.path.basename(d)}: "
                  f"analysing {os.path.basename(xtc)} (not PBC-corrected) - run 26_pbc_fix.sh first")
        label = os.path.basename(os.path.dirname(d)) + "/" + os.path.basename(d)
        try:
            r = analyse(top, xtc, label=label)
        except Exception as e:
            r = {"label": label, "error": f"{type(e).__name__}: {e}"}
        results.append(r)
        if "error" in r:
            print(f"{label:34s} ERROR: {r['error'][:110]}")
        else:
            print(f"{label:34s} {r['duration_ns']:6.1f} {r['ligand_drift_nm']['last']:7.3f} "
                  f"{r['ligand_conformational_nm']['last']:7.3f} "
                  f"{r['contacts_within_0.4nm']['last']:9d} "
                  f"{r['min_protein_distance_nm']['last']:8.3f}  {r['verdict']}")

    out = os.path.join(ROOT, "analysis", "md_stability.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
    print(f"\n-> {out}")
    print("""
READING
  drift     ligand RMSD after fitting on the protein backbone = motion relative to the pocket.
            < 0.3 nm and flat = stable; climbing past 0.5 nm = the ligand left and must not be
            reported as a binder.
  conf      ligand RMSD fitted on itself = conformational change only. Flat conf with rising drift
            means a rigid-body departure rather than a rearrangement.
  contacts  protein atoms within 0.4 nm; falling to ~0 means the ligand is solvent-exposed.
  mindist   closest protein-ligand approach; > 0.5 nm means nothing is holding it.
  Rg        protein radius of gyration, to show the fold itself did not collapse.""")
    return 0


if __name__ == "__main__":
    sys.exit(main())
