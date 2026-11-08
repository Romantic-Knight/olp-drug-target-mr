#!/usr/bin/env python3
"""fix_mte_bonds.py — repair the molybdopterin bond perception, then reparameterise.

The problem
-----------
tleap rejected MTE:  "Atom .R<MTE 1>.A<C4 5> has force field coordination 4 but only 3 bonded
neighbors."  Investigation of `mte.mol2` and of the crystal geometry found the cause is bond
perception, in three parts:

  1. **Three duplicate bonds** - antechamber wrote 1-12, 2-4 and 5-6 twice each. tleap sums bonds,
     so each duplicate steals a valence slot.
  2. **C4=O4 was written as a single bond.** Measured in 2A9A, C4-O4 is **1.249 A**, which is a
     carbonyl (C=O ~1.19-1.23 A), not a C-OH/C-OR single bond (~1.32-1.43 A). With it as a single
     bond and no hydrogen on C4, carbon 4 is left one bond short.
  3. Consequently the atom types are wrong too: C4 was typed `c3` (sp3) when it is sp2, and O4 was
     typed `os` (ether) when it is a carbonyl oxygen `o`.

This matters beyond tleap: the AM1-BCC charges were computed on that incorrect connectivity, so a
bond-order-only patch would leave the electrostatics wrong. The fix therefore rewrites the bond
orders and then **re-runs antechamber from a corrected SDF**, so types and charges are both derived
from the right structure.

The crystal geometry that justifies step 2 (from `2A9A.pdb`, chain A, MTE):
    C4-O4  1.249 A     -> C=O
    C4-N3  1.386 A     -> C-N
    C4-C9  1.396 A     -> C-C
three connections plus the carbonyl double bond = valence 4, which is what tleap expects.

Usage:
    python3 fix_mte_bonds.py <complex-dir> <net-charge>
"""
import os
import re
import subprocess
import sys

BOND_RE = re.compile(r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*$")


def fix_mol2(path):
    """Deduplicate bonds and promote the C4=O4 bond to double. Returns a short report."""
    lines = open(path, encoding="utf-8", errors="replace").read().split("\n")

    # map atom index -> atom name, so the bond is identified chemically rather than by number
    name_of, in_atoms = {}, False
    for ln in lines:
        if ln.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            continue
        if ln.startswith("@<TRIPOS>") and in_atoms:
            in_atoms = False
        elif in_atoms and ln.split():
            f = ln.split()
            name_of[int(f[0])] = f[1]

    out, seen, in_bonds = [], set(), False
    dupes, promoted = 0, None
    for ln in lines:
        if ln.startswith("@<TRIPOS>BOND"):
            in_bonds = True
            out.append(ln)
            continue
        if ln.startswith("@<TRIPOS>") and in_bonds:
            in_bonds = False

        if in_bonds:
            m = BOND_RE.match(ln)
            if m:
                i, j, order = int(m.group(2)), int(m.group(3)), int(m.group(4))
                key = (min(i, j), max(i, j))
                if key in seen:
                    dupes += 1
                    continue
                seen.add(key)
                pair = {name_of.get(i, "?"), name_of.get(j, "?")}
                if pair == {"C4", "O4"} and order == 1:
                    order = 2                      # the carbonyl, per the 1.249 A measurement
                    promoted = (i, j)
                    ln = f"{m.group(1):>6s}{i:>7d}{j:>7d}{order:>2d}\n"
        out.append(ln)

    open(path, "w", encoding="utf-8", newline="\n").write("\n".join(out))
    return dupes, promoted


def main():
    d = sys.argv[1]
    q = sys.argv[2] if len(sys.argv) > 2 else "-2"
    mol2 = os.path.join(d, "mte.mol2")
    if not os.path.exists(mol2):
        print(f"no {mol2}"); return 1

    print("=== 1. repair the mol2 ===")
    dupes, promoted = fix_mol2(mol2)
    print(f"  duplicate bonds removed : {dupes}")
    print(f"  C4=O4 promoted to order 2: {'yes' if promoted else 'NO - pair not found'}")

    # report the resulting valence of every atom, so over-coordination cannot slip through
    deg = {}
    for ln in open(mol2, encoding="utf-8", errors="replace"):
        m = BOND_RE.match(ln)
        if m:
            for i in (int(m.group(2)), int(m.group(3))):
                deg[i] = deg.get(i, 0) + 1
    over = sorted(i for i, n in deg.items() if n > 4)
    print(f"  atoms with more than 4 bonds: {over if over else 'none'}")

    print("\n=== 2. mol2 -> SDF so antechamber sees the corrected bond orders ===")
    sdf = os.path.join(d, "mte_fixed.sdf")
    r = subprocess.run(["obabel", mol2, "-O", sdf], capture_output=True, text=True)
    print(f"  obabel rc={r.returncode} {(r.stderr or '').strip()[:80]}")
    if not os.path.exists(sdf):
        print("  failed to write the SDF"); return 1

    print("\n=== 3. reparameterise from the corrected structure ===")
    for f in ("mte_fixed.mol2", "mte_fixed.frcmod"):
        p = os.path.join(d, f)
        if os.path.exists(p):
            os.remove(p)
    r = subprocess.run(["antechamber", "-i", "mte_fixed.sdf", "-fi", "sdf",
                        "-o", "mte_fixed.mol2", "-fo", "mol2", "-c", "bcc", "-s", "2",
                        "-at", "gaff2", "-nc", q, "-rn", "MTE"],
                       cwd=d, capture_output=True, text=True)
    print(f"  antechamber rc={r.returncode}")
    out_mol2 = os.path.join(d, "mte_fixed.mol2")
    if not os.path.exists(out_mol2):
        tail = open(os.path.join(d, "antechamber.log"), encoding="utf-8", errors="replace").read()[-500:]
        print("  FAILED; log tail:", tail)
        return 1

    r = subprocess.run(["parmchk2", "-i", "mte_fixed.mol2", "-f", "mol2",
                        "-o", "mte_fixed.frcmod", "-s", "gaff2"], cwd=d, capture_output=True, text=True)
    print(f"  parmchk2 rc={r.returncode}")

    # replace the broken files so the build script picks up the corrected ones
    for src, dst in (("mte_fixed.mol2", "mte.mol2"), ("mte_fixed.frcmod", "mte.frcmod")):
        s, t = os.path.join(d, src), os.path.join(d, dst)
        if os.path.exists(s):
            if os.path.exists(t):
                os.remove(t)
            os.rename(s, t)
    print("  mte.mol2 / mte.frcmod replaced with the corrected versions")

    # verification
    txt = open(os.path.join(d, "mte.mol2"), encoding="utf-8", errors="replace").read()
    atoms = 0
    in_atoms = False
    types = {}
    for ln in txt.split("\n"):
        if ln.startswith("@<TRIPOS>ATOM"):
            in_atoms = True; continue
        if ln.startswith("@<TRIPOS>") and in_atoms:
            in_atoms = False
        elif in_atoms and ln.split():
            f = ln.split()
            if len(f) >= 6:
                atoms += 1
                types[f[5]] = types.get(f[5], 0) + 1
    nb = sum(1 for ln in txt.split("\n") if BOND_RE.match(ln))
    print(f"\n  final: {atoms} atoms, {nb} bonds")
    print(f"  GAFF2 atom types in use: {dict(sorted(types.items()))}")
    frc = os.path.join(d, "mte.frcmod")
    if os.path.exists(frc):
        print(f"  frcmod revision flags: "
              f"{sum(1 for l in open(frc) if 'ATTN, need revision' in l)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
