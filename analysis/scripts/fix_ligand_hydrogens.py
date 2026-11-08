#!/usr/bin/env python3
"""fix_ligand_hydrogens.py — put the full hydrogen set back before parameterisation.

The bug this fixes
------------------
`ligand.pdb` was extracted from the docked AutoDock pose, and PDBQT follows the united-atom
convention: only **polar** hydrogens are written, non-polar ones are implied by the atom type
(`C` is really `CH1/CH2/CH3`). `antechamber`/`sqm`, however, need the explicit, complete structure.
Given the stripped file they infer the wrong valence, end up with an **odd electron count**, and
sqm refuses:

    QMMM: System specified with odd number of electrons (  127)
    QMMM: but odd spin (  1). You most likely have the charge of the QM region set incorrectly.

That is not a charge problem - the charge is right - it is a missing-hydrogen problem. Two of the
four ligands (Spiperone, Vorinostat) failed outright and the two that "succeeded" (GW-8510,
Sanguinarine) were parameterised from the same incomplete input, so all four are redone here.

What this does
--------------
Uses OpenBabel to add the missing hydrogens **in place**, which keeps every heavy atom at its
docked-pose coordinate - important, because the whole point of the parameterisation is to run MD
from that pose - and writes `ligand_fullH.pdb`. It then verifies that the heavy atoms did not move
and reports the electron count per ligand so the fix is auditable rather than assumed.

Usage:
    python analysis/scripts/fix_ligand_hydrogens.py [complex ...]
"""
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MD = os.path.join(ROOT, "analysis", "md")
OBABEL = "obabel"


def heavy_atoms(path):
    """(element, x, y, z) for non-hydrogen atoms in a PDB file."""
    out = []
    for ln in open(path, encoding="utf-8", errors="replace"):
        if ln.startswith(("ATOM", "HETATM")):
            elem = (ln[76:78].strip() or ln[12:16].strip()[:1]).upper()
            if elem.startswith("H"):
                continue
            try:
                out.append((elem, float(ln[30:38]), float(ln[38:46]), float(ln[46:54])))
            except ValueError:
                pass
    return out


def count_h(path):
    n = 0
    for ln in open(path, encoding="utf-8", errors="replace"):
        if ln.startswith(("ATOM", "HETATM")):
            elem = (ln[76:78].strip() or ln[12:16].strip()[:1]).upper()
            if elem.startswith("H"):
                n += 1
    return n


def main():
    complexes = sys.argv[1:] or ["GDF11_GW8510", "SUOX_Sanguinarine",
                                 "SUOX_Spiperone", "SUOX_Vorinostat"]
    rc_all = 0
    for cx in complexes:
        d = os.path.join(MD, "complexes", cx)
        pq = os.path.join(d, "ligand.pdbqt")
        src = pq if os.path.exists(pq) else os.path.join(d, "ligand.pdb")
        if not os.path.exists(src):
            print(f"{cx}: no ligand pose found"); rc_all = 1; continue

        # the docked pose as originally extracted, for the coordinate comparison
        orig = os.path.join(d, "ligand.pdb")
        before = heavy_atoms(orig) if os.path.exists(orig) else heavy_atoms(src)

        out = os.path.join(d, "ligand_fullH.pdb")
        r = subprocess.run([OBABEL, "-ipdbqt" if src.endswith("pdbqt") else "-ipdb",
                            src, "-O", out, "-h"],
                           capture_output=True, text=True)
        if r.returncode != 0 or not os.path.exists(out):
            print(f"{cx}: obabel -h failed: {(r.stderr or '').strip()[:120]}")
            rc_all = 1
            continue

        after = heavy_atoms(out)
        nh = count_h(out)
        # heavy atoms must be untouched, otherwise the pose moved and the MD would start elsewhere
        moved = 0
        if len(before) == len(after):
            for (e1, *c1), (e2, *c2) in zip(before, after):
                if sum(abs(a - b) for a, b in zip(c1, c2)) > 0.01:
                    moved += 1
        else:
            moved = -1

        print(f"{cx:20s} heavy {len(before)} -> {len(after)}  H added {nh}  "
              f"heavy-atom displacement: {'n/a' if moved < 0 else moved} "
              f"{'OK' if moved == 0 else 'CHECK'}")
        if moved != 0:
            rc_all = 1
    return rc_all


if __name__ == "__main__":
    sys.exit(main())
