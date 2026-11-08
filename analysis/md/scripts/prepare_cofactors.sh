#!/usr/bin/env bash
# prepare_cofactors.sh — GAFF2 parameters for the SUOX active-site cofactor and its sulfate.
#
# Why this script exists
# ----------------------
# The first pass built SUOX with a protein-only receptor, so the ligand sat in an emptied cavity.
# That was the single biggest caveat on the MD result. Two pieces of the real site were missing:
#
#   MTE  molybdopterin, the dithiolene cofactor that lines the pocket
#   SO4  the catalytic sulfate seen 5.1 A from the Mo in 2A9A
#
# MTE failed before for a specific reason: it was extracted from the PDB, which carries no bond
# orders, so antechamber had to perceive them - it produced three duplicate bonds and left the
# exocyclic C4-O4 as a single bond when the crystal distance is 1.249 A (a carbonyl). It then also
# inferred a net charge of +1 where the dithiolene/phosphate chemistry implies -2, and AM1-BCC
# charges built on the wrong charge state would misrepresent the electrostatics the ligand sees.
#
# Fixes applied here:
#   * feed an SDF with explicit bond orders (no perception needed -> no duplicate bonds)
#   * force the net charge with `-nc -2`. C10H16N5O6PS2 at -2 gives 208 electrons (even), which sqm
#     requires; the number is a *choice* and is recorded as such in the run notes.
#
# Mo is still absent, and deliberately: Vina cannot represent it, antechamber/mcpb cannot derive its
# parameters without a QM code, and MCPB.py accepts only Gaussian, which is not installed. Chemical
# accuracy here is therefore "cofactor present, metal absent", which is strictly better than an
# empty pocket but is not a complete metal site.
# NOTE: `set -u` must come AFTER sourcing the environment. activate.sh references variables that
# are unset, so with -u already active the script dies silently on the source line and prints
# nothing at all - the same failure that once made the PBC-correction step a no-op.
MD="/root/olp/analysis/md"
source /opt/md/activate.sh >/dev/null 2>&1
set -uo pipefail

FIXTURE="$MD/complexes/SUOX_Sanguinarine"
cd "$FIXTURE" || exit 1

echo "=== 1. MTE from the explicit-bond SDF, charge forced to -2 ==="
if [ ! -s mte_fullH.sdf ]; then
  echo "  FAIL mte_fullH.sdf missing or empty"; exit 1
fi
rm -f mte.mol2 mte.frcmod
antechamber -i mte_fullH.sdf -fi sdf -o mte.mol2 -fo mol2 \
            -c bcc -nc -2 -at gaff2 -rn MTE -pf y >mte_antechamber.log 2>&1
if [ -s mte.mol2 ]; then
  parmchk2 -i mte.mol2 -f mol2 -o mte.frcmod -s gaff2 >mte_parmchk2.log 2>&1
  echo "  mol2 $(stat -c%s mte.mol2) bytes, frcmod $(stat -c%s mte.frcmod) bytes"
else
  echo "  FAIL - see mte_antechamber.log"; tail -6 mte_antechamber.log
fi

echo
echo "=== 2. catalytic sulfate (SO4 7503, 5.1 A from Mo) from 2A9A ==="
python3 - <<'PY'
import os, math
src = "/root/olp/analysis/md/complexes/SUOX_Sanguinarine/../2A9A.pdb"
# the trimmed receptor is what the MD uses; pull the sulfate out of the original crystal file
cand = ["/root/olp/analysis/docking/structures/2A9A.pdb",
        "/root/olp/analysis/md/structures/2A9A.pdb", "2A9A.pdb", "../2A9A.pdb"]
p = next((c for c in cand if os.path.exists(c)), None)
if p is None:
    raise SystemExit("2A9A.pdb not found; cannot extract SO4")
mo, sulfates = None, {}
for l in open(p):
    if not l.startswith(("ATOM", "HETATM")):
        continue
    if l[21] != "A":
        continue
    res = l[17:20].strip(); rid = l[22:26].strip()
    xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    if res == "MO":
        mo = xyz
    elif res == "SO4":
        sulfates.setdefault(rid, []).append((l[12:16], xyz))
if mo is None or not sulfates:
    raise SystemExit(f"MO found={mo is not None}, sulfates={list(sulfates)}")
# keep the sulfate nearest the metal: that is the catalytic one
best = min(sulfates, key=lambda k: min(math.dist(mo, a[1]) for a in sulfates[k]))
out = ["REMARK catalytic sulfate from 2A9A chain A, 5.1 A from the Mo",
       "HETATM"]
# renumber and rename to SO4 for tleap
lines = []
for i, (an, xyz) in enumerate(sulfates[best], start=1):
    lines.append(f"HETATM{i:5d} {an:>4s} SO4  A 751   " +
                 f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00          S\n")
open("so4.pdb", "w").write("\n".join(lines) + "END\n")
print(f"  wrote so4.pdb: {len(lines)} atoms (resid {best}, "
      f"{min(math.dist(mo, a[1]) for a in sulfates[best]):.2f} A from Mo)")
PY

if [ -s so4.pdb ]; then
  obabel so4.pdb -O so4.sdf --gen3d >/dev/null 2>&1 || obabel so4.pdb -O so4.sdf >/dev/null 2>&1
  rm -f so4.mol2 so4.frcmod
  antechamber -i so4.sdf -fi sdf -o so4.mol2 -fo mol2 -c bcc -nc -2 -at gaff2 -rn SO4 -pf y >so4_antechamber.log 2>&1
  if [ -s so4.mol2 ]; then
    parmchk2 -i so4.mol2 -f mol2 -o so4.frcmod -s gaff2 >so4_parmchk2.log 2>&1
    echo "  mol2 $(stat -c%s so4.mol2) bytes, frcmod $(stat -c%s so4.frcmod) bytes"
  else
    echo "  SO4 FAIL - see so4_antechamber.log"; tail -6 so4_antechamber.log
  fi
fi

echo
echo "=== 3. verification ==="
python3 - <<'PY'
import os, collections
def check(name):
    p = f"{name}.mol2"
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        print(f"  {name}: MOL2 MISSING"); return
    sec, atoms, bonds = None, {}, []
    for l in open(p):
        if l.startswith("@<TRIPOS>"):
            sec = l.strip(); continue
        if sec == "@<TRIPOS>ATOM":
            f = l.split()
            if len(f) >= 6:
                atoms[int(f[0])] = (f[1], f[5])
        elif sec == "@<TRIPOS>BOND":
            f = l.split()
            if len(f) >= 4:
                bonds.append((int(f[1]), int(f[2]), f[3]))
    try:
        chg = [float(v[1]) for v in atoms.values()]
    except ValueError:
        chg = []
    pairs = [tuple(sorted((a, b))) for a, b, _ in bonds]
    dup = [k for k, c in collections.Counter(pairs).items() if c > 1]
    print(f"  {name}: atoms={len(atoms)} bonds={len(bonds)} "
          f"charge={sum(chg):+.2f} duplicate_bonds={len(dup)}")
    if dup:
        print(f"      duplicated: {dup[:5]}")
for n in ("mte", "so4"):
    check(n)
PY
echo
echo "=== 4. what tleap needs to know ==="
echo "  loadamberparams mte.frcmod / so4.frcmod"
echo "  MTE = loadmol2 mte.mol2 ; SO4 = loadmol2 so4.mol2"
