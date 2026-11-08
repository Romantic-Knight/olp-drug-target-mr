#!/usr/bin/env bash
# 32_build_gdf11.sh — build and equilibrate the GDF11 (5E4G) monomer system.
#
# Three things have to be handled before this receptor can go into Amber, and none of them were
# problems for SUOX:
#
# 1. **Alternate locations.** Residues 5, 11, 16, 73, 75 and 90 carry A and B conformers (42
#    duplicate atom records). tleap cannot take altLoc characters and would either reject the
#    residue or build it twice, so only the blank/A records are kept. A is the higher-occupancy
#    conformer in every case (0.73/0.27 and 0.50/0.50).
#
# 2. **Disulfides.** GDF11 is a TGF-beta-family ligand with a cysteine-knot fold: SG-SG pairs are
#    present at 2.0 A between Cys6-Cys16, Cys15-Cys74, Cys43-Cys106 and Cys47-Cys108. tleap only
#    forms a disulfide if it is asked to, so these four are written as explicit `bond` commands.
#    Without them every cysteine becomes a reduced thiol, which would be a real chemical error in
#    a fold held together by those bonds.
#
# 3. **Cys73 is the inter-chain cysteine.** In the crystal it is a symmetry-related disulfide to the
#    partner monomer (SG-SG 2.4 A within one asymmetric unit). This build is the monomer as agreed,
#    so Cys73 is left as a free thiol and that is recorded as a limitation.
#
# The ligand has no metal and no cofactor, so this uses the plain (non-holo) route.
#
# Usage: bash 32_build_gdf11.sh
MD="/root/olp/analysis/md"
source /opt/md/activate.sh >/dev/null 2>&1
set -uo pipefail
cd "$MD"

cx="${1:-GDF11_GW8510}"
dir="$MD/complexes/$cx"
MARGIN="${MARGIN:-1.2}"

echo "############ $cx"
[ -s "$dir/ligand.mol2" ] || { echo "  no ligand.mol2"; exit 1; }
[ -s "$dir/ligand.frcmod" ] || { echo "  no ligand.frcmod"; exit 1; }

echo
echo "=== 1. strip alternate locations, keep blank/A ==="
python3 - "$dir" <<'PY'
import os, sys, collections
d = sys.argv[1]
src = os.path.join(d, "receptor.pdb")
kept, dropped = [], collections.Counter()
for l in open(src):
    if l.startswith(("ATOM", "HETATM")):
        alt = l[16]
        if alt in (" ", "A"):
            # blank the altLoc column so no downstream tool sees it
            kept.append(l[:16] + " " + l[17:].rstrip("\n"))
        else:
            dropped[l[17:20].strip() + l[22:27].strip()] += 1
    elif l.startswith("TER") or l.startswith("END"):
        kept.append(l.rstrip("\n"))
out = os.path.join(d, "receptor_clean.pdb")
open(out, "w").write("\n".join(kept) + "\n")
print(f"  kept {len([k for k in kept if k.startswith('ATOM')])} atoms, "
      f"dropped {sum(dropped.values())} alternate-conformer records")
if dropped:
    print("  dropped from residues:", ", ".join(sorted(dropped)))
PY

echo
echo "=== 2. find disulfides on the cleaned structure, write tleap bonds ==="
python3 - "$dir" <<'PY'
import os, sys, math
d = sys.argv[1]
sg = []
for l in open(os.path.join(d, "receptor_clean.pdb")):
    if l.startswith("ATOM") and l[12:16].strip() == "SG" and l[17:20].strip() == "CYS":
        sg.append((l[22:26].strip(), (float(l[30:38]), float(l[38:46]), float(l[46:54]))))
lines, used = [], set()
for i in range(len(sg)):
    for j in range(i + 1, len(sg)):
        a, b = sg[i], sg[j]
        dist = math.dist(a[1], b[1])
        if dist < 2.5 and a[0] not in used and b[0] not in used:
            lines.append((a[0], b[0], dist))
            used.update((a[0], b[0]))
with open(os.path.join(d, "disulfides.txt"), "w") as f:
    for a, b, dist in lines:
        f.write(f"bond rec.{int(a)}.SG rec.{int(b)}.SG\n")
    f.write(f"# {len(lines)} disulfides; Cys73 is the inter-chain cysteine and is left reduced\n")
print(f"  {len(lines)} disulfides:")
for a, b, dist in lines:
    print(f"    Cys{a} - Cys{b}   {dist:.2f} A")
free = [r for r, _ in sg if r not in used]
print(f"  free thiols (no partner in this monomer): {', '.join(free) if free else 'none'}")

# Amber convention: a cysteine in a disulfide is residue CYX, whose SG is atom type S (thioether).
# Left as CYS its SG is typed SH (thiol) and tleap reports
#   "Could not find bond parameter for atom types: SH - SH"
# so every paired cysteine is renamed here and the unpaired one stays CYS.
cyx = {a for a, _ in sg if a in used} | {b for _, b, _ in lines}
cyx |= {r for r, _ in sg if r not in free}
src = os.path.join(d, "receptor_clean.pdb")
out = []
changed = set()
for l in open(src):
    if l.startswith("ATOM") and l[17:20] == "CYS" and l[22:26].strip() in cyx:
        changed.add(l[22:26].strip())
        l = l[:17] + "CYX" + l[20:]
    out.append(l.rstrip("\n"))
open(src, "w").write("\n".join(out) + "\n")
print(f"  renamed to CYX: {', '.join(sorted(changed, key=int))}")
PY

echo
echo "=== 3. tleap with explicit disulfide bonds ==="
cat > "$dir/leap_gdf11.in" <<EOF
source leaprc.protein.ff19SB
source leaprc.gaff2
source leaprc.water.tip3p
loadamberparams ligand.frcmod
LIG = loadmol2 ligand.mol2
rec = loadpdb receptor_clean.pdb
$(cat "$dir/disulfides.txt" | grep -v '^#')
cx = combine { rec LIG }
check cx
saveamberparm cx gdf11_dry.prmtop gdf11_dry.inpcrd
solvateoct cx TIP3PBOX $MARGIN
addions cx Na+ 0
addions cx Cl- 0
saveamberparm cx gdf11_solv.prmtop gdf11_solv.inpcrd
savepdb cx gdf11_solv.pdb
quit
EOF
( cd "$dir" && tleap -f leap_gdf11.in >tleap_gdf11.log 2>&1 )
if [ ! -s "$dir/gdf11_solv.prmtop" ]; then
    echo "  tleap FAILED"; grep -iE 'fatal|error|does not have' "$dir/tleap_gdf11.log" | head -6 | sed 's/^/      /'
    exit 1
fi
dry=$(awk '/%FLAG POINTERS/{getline; getline; print $1}' "$dir/gdf11_dry.prmtop")
sol=$(awk '/%FLAG POINTERS/{getline; getline; print $1}' "$dir/gdf11_solv.prmtop")
echo "  tleap OK: dry ${dry} atoms, solvated ${sol} atoms"

echo
echo "=== 4. parmed -> GROMACS ==="
cat > "$dir/to_gmx_gdf11.py" <<'PY'
from parmed.amber import AmberParm
from parmed.gromacs import GromacsTopologyFile
p = AmberParm("gdf11_solv.prmtop", "gdf11_solv.inpcrd")
g = GromacsTopologyFile.from_structure(p)
g.save("gdf11.top", overwrite=True)
p.save("gdf11.gro", overwrite=True)
print("converted", len(p.atoms), "atoms")
PY
( cd "$dir" && python to_gmx_gdf11.py >to_gmx_gdf11.log 2>&1 )
[ -s "$dir/gdf11.top" ] || { echo "  parmed FAILED"; tail -4 "$dir/to_gmx_gdf11.log" | sed 's/^/      /'; exit 1; }
echo "  gdf11.top written"

echo
echo "=== 5. EM -> NVT -> NPT ==="
for m in em nvt npt; do cp -n "$dir/../SUOX_Sanguinarine/$m.mdp" "$dir/$m.mdp" 2>/dev/null || true; done
( cd "$dir"
  gmx grompp -f em.mdp  -c gdf11.gro -p gdf11.top -o gem.tpr -maxwarn 2 >gg_em.log 2>&1 || exit 11
  gmx mdrun -deffnm gem -nb gpu >mg_em.log 2>&1 || gmx mdrun -deffnm gem >mg_em.log 2>&1 || exit 12
  gmx grompp -f nvt.mdp -c gem.gro -r gem.gro -p gdf11.top -o gnvt.tpr -maxwarn 2 >gg_nvt.log 2>&1 || exit 13
  gmx mdrun -deffnm gnvt -nb gpu -pme gpu -bonded gpu -update gpu -ntomp 4 >mg_nvt.log 2>&1 || \
    gmx mdrun -deffnm gnvt >mg_nvt.log 2>&1 || exit 14
  gmx grompp -f npt.mdp -c gnvt.gro -r gnvt.gro -t gnvt.cpt -p gdf11.top -o gnpt.tpr -maxwarn 2 >gg_npt.log 2>&1 || exit 15
  gmx mdrun -deffnm gnpt -nb gpu -pme gpu -bonded gpu -update gpu -ntomp 4 >mg_npt.log 2>&1 || \
    gmx mdrun -deffnm gnpt >mg_npt.log 2>&1 || exit 16 )
rc=$?
if [ "$rc" -ne 0 ]; then
    echo "  equilibration FAILED at step $rc"
    for l in gg_em gg_nvt gg_npt mg_em mg_nvt mg_npt; do
        [ -f "$dir/$l.log" ] && grep -iE 'error|fatal' "$dir/$l.log" | head -2 | sed "s/^/      $l: /"
    done
    exit 1
fi
echo "  equilibrated (EM, NVT, NPT)"
echo "$cx READY_GDF11 dry=$dry solv=$sol" >> "$MD/PRODUCTION_GDF11_STATUS.txt"
cat "$MD/PRODUCTION_GDF11_STATUS.txt"
