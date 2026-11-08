#!/usr/bin/env python3
"""interactions_md.py — hydrogen bonds and residue-level contacts from the harvested trajectories.

Produces the numbers behind the two interaction panels:
  * per (complex, residue) contact frequency  -> hotspot map
  * per (complex, residue) hydrogen-bond occupancy -> H-bond heat map
  * hydrogen-bond count per frame -> distribution per complex

Geometry follows the standard GROMACS/PLIP criteria:
  hydrogen bond   donor-H...acceptor with H...A <= 0.25 nm and angle(D-H...A) >= 130 deg
  contact         any protein heavy atom within 0.4 nm of any ligand heavy atom

MDAnalysis reports Angstrom, so every threshold below is multiplied by 10 explicitly rather than
being written as a nm value (that mistake produced 10x-wrong contact counts earlier in this project).

Outputs
  results/tables/md_residue_contacts.csv
  results/tables/md_hbonds.csv
"""
import glob
import os
import sys
from collections import defaultdict

import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance, calc_angles

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
BASE = os.path.join(ROOT, "remote", "easyai", "md_solute_20261018", "unpacked")

CONTACT_A = 4.0          # 0.4 nm
HB_DIST_A = 2.5          # 0.25 nm  H...acceptor
HB_ANGLE_DEG = 130.0
HB_ANGLE_RAD = np.radians(130.0)   # calc_angles returns radians, not degrees
STRIDE = 10              # every 1 ns for a 100 ps trajectory


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
    lig = u.select_atoms("resname LIG or resname UNL")
    if len(lig) == 0:
        return None
    prot = u.select_atoms("protein")
    # polar atoms only for the H-bond search
    donors = u.select_atoms("protein and (name N* or name O* or name S*)")
    lig_n = lig.select_atoms("name N* or name O* or name S*")

    nframes = 0
    res_contact = defaultdict(int)     # resid -> frames with >=1 contact
    res_hb = defaultdict(int)          # resid -> frames with >=1 hbond
    hb_counts = []

    # donor hydrogens: GROMACS/Amber names start with H, bonded to N/O/S
    hyd = u.select_atoms("protein and name H*")
    for ts in u.trajectory[::STRIDE]:
        nframes += 1
        # --- contacts: residue level ---
        pairs, dists = capped_distance(prot.positions, lig.positions,
                                       max_cutoff=CONTACT_A, return_distances=True)
        if len(pairs):
            hit = np.unique(prot.resindices[pairs[:, 0]])
            for ri in hit:
                res_contact[u.residues[ri].resid] += 1

        # --- hydrogen bonds, both directions ---
        # The ligand here carries only two polar hydrogens (H7 on N1, H16 on N2 for Spiperone), so
        # searching protein-donor -> ligand-acceptor alone reports zero. A protein donor -> ligand
        # acceptor AND a ligand donor -> protein acceptor search are both required.
        # A .gro has no bond table, so each donor is found as the nearest covalently bound N/O/S.
        nhb = 0

        def _donor_heavy(h_pos, heavy_pos, heavy_names):
            d = np.linalg.norm(heavy_pos - h_pos, axis=1)
            j = int(np.argmin(d))
            if d[j] < 1.3 and heavy_names[j][0] in ("N", "O", "S"):
                return j
            return None

        prot_heavy = prot.select_atoms("name N* or name O* or name S*")
        pnames = list(prot_heavy.names)

        # (i) protein donates, ligand accepts
        if len(hyd) and len(lig_n):
            hpairs, _ = capped_distance(hyd.positions, lig_n.positions,
                                        max_cutoff=HB_DIST_A, return_distances=True)
            for hi, ai in hpairs:
                h_atom = u.atoms[hyd.indices[hi]]
                j = _donor_heavy(h_atom.position, prot_heavy.positions, pnames)
                if j is None:
                    continue
                acc = u.atoms[lig_n.indices[ai]]
                ang = calc_angles(prot_heavy.positions[j][None, :], h_atom.position[None, :],
                                  acc.position[None, :], backend="serial")[0]
                if ang >= HB_ANGLE_RAD:
                    nhb += 1
                    res_hb[h_atom.residue.resid] += 1

        # (ii) ligand donates, protein accepts
        lig_h = lig.select_atoms("name H*")
        lig_heavy = lig.select_atoms("not name H*")
        lnames = list(lig_heavy.names)
        if len(lig_h) and len(prot_heavy):
            hpairs, _ = capped_distance(lig_h.positions, prot_heavy.positions,
                                        max_cutoff=HB_DIST_A, return_distances=True)
            for hi, ai in hpairs:
                h_atom = u.atoms[lig_h.indices[hi]]
                j = _donor_heavy(h_atom.position, lig_heavy.positions, lnames)
                if j is None:
                    continue
                acc = u.atoms[prot_heavy.indices[ai]]
                ang = calc_angles(lig_heavy.positions[j][None, :], h_atom.position[None, :],
                                  acc.position[None, :], backend="serial")[0]
                if ang >= HB_ANGLE_RAD:
                    nhb += 1
                    res_hb[acc.residue.resid] += 1
        hb_counts.append(nhb)

    return dict(nframes=nframes,
                res_contact={k: v / nframes for k, v in res_contact.items()},
                res_hb={k: v / nframes for k, v in res_hb.items()},
                hb_mean=float(np.mean(hb_counts)) if hb_counts else 0.0,
                hb_max=int(np.max(hb_counts)) if hb_counts else 0)


def main():
    rs = runs()
    print(f"runs: {len(rs)}")
    os.makedirs(os.path.join(ROOT, "results", "tables"), exist_ok=True)

    contact_rows, hb_rows = [], []
    for cx, run, top, xtc in rs:
        r = analyse(top, xtc)
        if r is None:
            continue
        label = f"{cx}/{run}"
        for resid, f in sorted(r["res_contact"].items()):
            contact_rows.append(dict(run=label, resid=resid, contact_freq=round(f, 3)))
        for resid, f in sorted(r["res_hb"].items()):
            hb_rows.append(dict(run=label, resid=resid, hbond_occ=round(f, 3)))
        print(f"  {label:44s} frames {r['nframes']:4d}  "
              f"residues contacted {len(r['res_contact']):3d}  "
              f"H-bond mean {r['hb_mean']:.1f} max {r['hb_max']}")

    import csv
    p1 = os.path.join(ROOT, "results", "tables", "md_residue_contacts.csv")
    with open(p1, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["run", "resid", "contact_freq"])
        w.writeheader(); w.writerows(contact_rows)
    p2 = os.path.join(ROOT, "results", "tables", "md_hbonds.csv")
    with open(p2, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["run", "resid", "hbond_occ"])
        w.writeheader(); w.writerows(hb_rows)
    print(f"\nwrote {os.path.relpath(p1, ROOT)} ({len(contact_rows)} rows)")
    print(f"wrote {os.path.relpath(p2, ROOT)} ({len(hb_rows)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
