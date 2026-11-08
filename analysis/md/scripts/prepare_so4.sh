#!/usr/bin/env bash
# prepare_so4.sh — sulfate parameters, taking the direct PDB route.
#
# Two bugs are being worked around here:
#
# 1. The first attempt ran `obabel so4.pdb -O so4.sdf --gen3d` and then antechamber on the SDF. It
#    failed with "Invalid number of atoms (-999); MDL SDF supports at most 999 atoms", i.e. the SDF
#    was malformed. The cause was the hand-written PDB: the element symbol was appended after the
#    temperature-factor field instead of into columns 77-78, so nothing parsed it as an element.
#
# 2. antechamber reads PDB directly (`-fi pdb`), which skips the conversion step that introduced the
#    malformed intermediate. Bond orders for SO4 are unambiguous, so no SDF is needed at all.
#
# SO4 carries -2: S(VI) with four oxygens and no protons.
MD="/root/olp/analysis/md"
source /opt/md/activate.sh >/dev/null 2>&1
set -uo pipefail

FIXTURE="$MD/complexes/SUOX_Sanguinarine"
cd "$FIXTURE" || exit 1

echo "=== rewrite so4.pdb with correct PDB column layout ==="
python3 - <<'PY'
import os, math
cand = ["/root/olp/analysis/docking/structures/2A9A.pdb",
        "/root/olp/analysis/md/structures/2A9A.pdb"]
p = next((c for c in cand if os.path.exists(c)), None)
if p is None:
    raise SystemExit("2A9A.pdb not found")

mo, sulfates = None, {}
for l in open(p):
    if not l.startswith(("ATOM", "HETATM")) or l[21] != "A":
        continue
    res, rid = l[17:20].strip(), l[22:26].strip()
    xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    if res == "MO":
        mo = xyz
    elif res == "SO4":
        sulfates.setdefault(rid, [(l[12:16].strip(), xyz)])
        sulfates[rid].append((l[12:16].strip(), xyz))
# keep only the sulfate nearest the metal
best = min(sulfates, key=lambda k: min(math.dist(mo, a[1]) for a in sulfates[k]))
atoms = sulfates[best][1:] if len(sulfates[best]) > 5 else sulfates[best]

out = ["REMARK  catalytic sulfate from 2A9A chain A (nearest the Mo)"]
for i, (an, (x, y, z)) in enumerate(atoms, start=1):
    elem = an[0] if an[0].isalpha() else an[1]
    # columns: 1-6 record, 7-11 serial, 13-16 name, 18-20 resName, 22 chain,
    #          23-26 resSeq, 31-38 x, 39-46 y, 47-54 z, 55-60 occ, 61-66 bfac, 77-78 element
    line = ("HETATM" + f"{i:5d}" + " " + f"{an:<4s}" + " " + "SO4" + " " + "A" +
            f"{751:4d}" + "    " + f"{x:8.3f}{y:8.3f}{z:8.3f}" +
            f"{1.00:6.2f}{0.00:6.2f}" + " " * 10 + f"{elem:>2s}")
    out.append(line)
out.append("END")
open("so4.pdb", "w").write("\n".join(out) + "\n")
print(f"  {len(atoms)} atoms, resSeq 751, nearest-atom distance to Mo "
      f"{min(math.dist(mo, a[1]) for a in atoms):.2f} A")
PY
echo "  --- file ---"
cat -A so4.pdb | head -3 | cut -c1-90

echo
echo "=== antechamber straight from the PDB ==="
rm -f so4.mol2 so4.frcmod
antechamber -i so4.pdb -fi pdb -o so4.mol2 -fo mol2 -c bcc -nc -2 -at gaff2 -rn SO4 -pf y >so4_antechamber.log 2>&1
if [ -s so4.mol2 ]; then
  parmchk2 -i so4.mol2 -f mol2 -o so4.frcmod -s gaff2 >so4_parmchk2.log 2>&1
  echo "  mol2 $(stat -c%s so4.mol2) bytes, frcmod $(stat -c%s so4.frcmod) bytes"
else
  echo "  FAIL"; tail -10 so4_antechamber.log
fi

echo
echo "=== verification (charge is field 9 of a mol2 ATOM line, not the atom type) ==="
python3 - <<'PY'
import os, collections
for name in ("mte", "so4"):
    p = f"{name}.mol2"
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        print(f"  {name:4s} MISSING"); continue
    sec, atoms, bonds = None, [], []
    for l in open(p):
        if l.startswith("@<TRIPOS>"):
            sec = l.strip(); continue
        f = l.split()
        if sec == "@<TRIPOS>ATOM" and len(f) >= 9:
            atoms.append((f[0], f[1], f[5], float(f[8])))
        elif sec == "@<TRIPOS>BOND" and len(f) >= 4:
            bonds.append((int(f[1]), int(f[2]), f[3]))
    chg = sum(a[3] for a in atoms)
    pairs = [tuple(sorted((a, b))) for a, b, _ in bonds]
    dup = [k for k, c in collections.Counter(pairs).items() if c > 1]
    types = collections.Counter(a[2] for a in atoms)
    print(f"  {name:4s} atoms={len(atoms):3d} bonds={len(bonds):3d} "
          f"charge={chg:+.2f} duplicate_bonds={len(dup)}")
    print(f"       types: {dict(types)}")
PY
