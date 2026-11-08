#!/usr/bin/env bash
# 27_build_suox_holo.sh — build SUOX with the cofactor present (MTE + catalytic sulfate).
#
# What changed from the first pass
# --------------------------------
# 23_build_suox_plain.sh deliberately dropped the cofactor because MTE could not be parameterised
# correctly (see its header). That is now fixed: MTE comes from an SDF with explicit bond orders
# and a forced net charge of -2, giving 40 atoms / 42 bonds / -2.00 / zero duplicate bonds, and the
# catalytic sulfate is parameterised the same way. So the ligand now sits in a pocket that contains
# the cofactor rather than an empty cavity - which was the largest caveat on the earlier result.
#
# Still absent, and unavoidably so: **Mo**. Vina cannot represent it, antechamber cannot type it,
# and MCPB.py accepts only Gaussian, which is not installed here. "Cofactor present, metal absent"
# is closer to the holoenzyme than before but is not a complete metal site.
#
# Assembly. The three pieces must share one coordinate frame, so they are written into a single PDB
# (protein from the trimmed receptor, MTE and SO4 straight from their mol2 coordinates, which came
# from the crystal) and loaded with loadpdb; the ligand is combined in separately, as in the plain
# build. Atom names in the PDB are taken verbatim from the mol2 files so that loadpdb can match the
# MTE and SO4 units loaded by loadmol2.
#
# Usage:  source /opt/md/activate.sh && bash 27_build_suox_holo.sh [complex ...]
MD="/root/olp/analysis/md"
source /opt/md/activate.sh >/dev/null 2>&1
set -uo pipefail
cd "$MD"

COMPLEXES="${*:-SUOX_Sanguinarine SUOX_Spiperone SUOX_Vorinostat}"
MARGIN="${MARGIN:-1.2}"

: > "$MD/PRODUCTION_HOLO_STATUS.txt"

for cx in $COMPLEXES; do
    dir="$MD/complexes/$cx"
    [ -d "$dir" ] || { echo "skip $cx (no dir)"; continue; }
    [ -s "$dir/ligand.mol2" ] || { echo "skip $cx (ligand not parameterised)"; continue; }
    echo "############ $cx"

    # the cofactor is identical in all three complexes (same crystal, same MTE/SO4), so the
    # parameter files are copied rather than re-derived. This must happen BEFORE the existence
    # check below - an earlier ordering made the check fire first, so the other two complexes were
    # skipped even though the files were available.
    for f in mte.mol2 mte.frcmod so4.mol2 so4.frcmod; do
        if [ ! -s "$dir/$f" ] && [ -s "$MD/complexes/SUOX_Sanguinarine/$f" ]; then
            cp "$MD/complexes/SUOX_Sanguinarine/$f" "$dir/$f"
            echo "  copied $f from SUOX_Sanguinarine"
        fi
    done
    if [ ! -s "$dir/mte.mol2" ] || [ ! -s "$dir/so4.mol2" ]; then
        echo "  skip: cofactor parameters still missing"
        continue
    fi

    # ---- assemble the holo receptor PDB -------------------------------------------------------
    python3 - "$dir" <<'PY'
import os, sys
d = sys.argv[1]

def mol2_atoms(path):
    """(name, x, y, z, element) in file order, so PDB and mol2 agree exactly."""
    out, sec = [], None
    for l in open(path):
        if l.startswith("@<TRIPOS>"):
            sec = l.strip(); continue
        if sec == "@<TRIPOS>ATOM":
            f = l.split()
            if len(f) >= 6:
                nm = f[1]
                elem = nm[0] if nm[0].isalpha() else (nm[1] if len(nm) > 1 else "C")
                out.append((nm, float(f[2]), float(f[3]), float(f[4]), elem))
    return out

prot = [l for l in open(os.path.join(d, "receptor_protein_only.pdb")) if l.startswith("ATOM")]
lines = ["REMARK  holo receptor: protein + MTE (molybdopterin) + catalytic SO4, no Mo",
         "REMARK  MTE and SO4 atom names and coordinates are taken verbatim from their mol2 files"]
lines += [l.rstrip("\n") for l in prot]
# TER closes the protein chain. Without it the cofactor HETATM records continue the same chain and
# tleap stops treating the last protein residue as C-terminal, which leaves its OXT untyped
# ("Atom .R<ASP 466>.A<OXT 13> does not have a type").
lines.append("TER")

# MTE and SO4 carry their crystal coordinates, so every piece stays in one frame
for resname, mol2, resseq in (("MTE", "mte.mol2", 7501), ("SO4", "so4.mol2", 751)):
    atoms = mol2_atoms(os.path.join(d, mol2))
    for i, (nm, x, y, z, elem) in enumerate(atoms, start=1):
        # PDB columns: 1-6 HETATM, 7-11 serial, 13-16 name, 18-20 resName, 22 chain,
        # 23-26 resSeq, 31-38 x, 39-46 y, 47-54 z, 55-60 occ, 61-66 bfac, 77-78 element
        lines.append("HETATM" + f"{i:5d}" + " " + f"{nm:<4s}" + " " + f"{resname:<3s}" + " " + "A" +
                     f"{resseq:4d}" + "    " + f"{x:8.3f}{y:8.3f}{z:8.3f}" +
                     f"{1.00:6.2f}{0.00:6.2f}" + " " * 10 + f"{elem:>2s}")
    print(f"  {resname}: {len(atoms)} atoms added at resSeq {resseq}")
lines.append("END")
open(os.path.join(d, "receptor_holo.pdb"), "w").write("\n".join(lines) + "\n")
cof = sum(1 for l in lines if l.startswith("HETATM"))
print(f"  receptor_holo.pdb: {len(prot)} protein + {cof} cofactor atoms, TER after protein")
print("  last protein atoms:", " | ".join(x[12:27].strip() for x in [l.rstrip() for l in prot[-3:]]))
PY

    # ---- tleap -------------------------------------------------------------------------------
    cat > "$dir/leap_holo.in" <<EOF
source leaprc.protein.ff19SB
source leaprc.gaff2
source leaprc.water.tip3p
loadamberparams ligand.frcmod
loadamberparams mte.frcmod
loadamberparams so4.frcmod
MTE = loadmol2 mte.mol2
SO4 = loadmol2 so4.mol2
LIG = loadmol2 ligand.mol2
rec = loadpdb receptor_holo.pdb
cx  = combine { rec LIG }
check cx
saveamberparm cx holo_dry.prmtop holo_dry.inpcrd
solvateoct cx TIP3PBOX $MARGIN
addions cx Na+ 0
addions cx Cl- 0
saveamberparm cx holo_solv.prmtop holo_solv.inpcrd
savepdb cx holo_solv.pdb
quit
EOF
    ( cd "$dir" && tleap -f leap_holo.in >tleap_holo.log 2>&1 )
    if [ ! -s "$dir/holo_solv.prmtop" ]; then
        echo "  tleap FAILED"
        grep -iE 'fatal|error|could not|does not have' "$dir/tleap_holo.log" | head -6 | sed 's/^/      /'
        continue
    fi
    atoms=$(awk '/%FLAG POINTERS/{getline; getline; print $1}' "$dir/holo_solv.prmtop")
    dry=$(awk '/%FLAG POINTERS/{getline; getline; print $1}' "$dir/holo_dry.prmtop")
    echo "  tleap OK: dry ${dry} atoms, solvated ${atoms} atoms"

    # ---- parmed -> GROMACS -------------------------------------------------------------------
    cat > "$dir/to_gmx_holo.py" <<'PY'
from parmed.amber import AmberParm
from parmed.gromacs import GromacsTopologyFile
p = AmberParm("holo_solv.prmtop", "holo_solv.inpcrd")
g = GromacsTopologyFile.from_structure(p)
g.save("holo.top", overwrite=True)
p.save("holo.gro", overwrite=True)
print("converted", len(p.atoms), "atoms")
PY
    ( cd "$dir" && python to_gmx_holo.py >to_gmx_holo.log 2>&1 )
    if [ ! -s "$dir/holo.top" ]; then
        echo "  parmed conversion FAILED"; tail -4 "$dir/to_gmx_holo.log" | sed 's/^/      /'; continue
    fi
    echo "  GROMACS files written (holo.top, holo.gro)"

    # ---- reuse the plain build's mdp files, or write them if absent ---------------------------
    for m in em nvt npt; do
        [ -s "$dir/$m.mdp" ] || { echo "  missing $m.mdp (run the plain build once first)"; break; }
    done
    if [ ! -s "$dir/em.mdp" ]; then continue; fi

    # ---- EM -> NVT -> NPT -------------------------------------------------------------------
    ( cd "$dir"
      gmx grompp -f em.mdp  -c holo.gro -p holo.top -o hem.tpr -maxwarn 2 >gh_em.log 2>&1 || exit 11
      gmx mdrun -deffnm hem -nb gpu >mh_em.log 2>&1 || gmx mdrun -deffnm hem >mh_em.log 2>&1 || exit 12
      gmx grompp -f nvt.mdp -c hem.gro -r hem.gro -p holo.top -o hnvt.tpr -maxwarn 2 >gh_nvt.log 2>&1 || exit 13
      gmx mdrun -deffnm hnvt -nb gpu -pme gpu -bonded gpu -update gpu >mh_nvt.log 2>&1 || \
        gmx mdrun -deffnm hnvt >mh_nvt.log 2>&1 || exit 14
      gmx grompp -f npt.mdp -c hnvt.gro -r hnvt.gro -t hnvt.cpt -p holo.top -o hnpt.tpr -maxwarn 2 >gh_npt.log 2>&1 || exit 15
      gmx mdrun -deffnm hnpt -nb gpu -pme gpu -bonded gpu -update gpu >mh_npt.log 2>&1 || \
        gmx mdrun -deffnm hnpt >mh_npt.log 2>&1 || exit 16 )
    rc=$?
    if [ $rc -ne 0 ]; then
        echo "  equilibration FAILED at step $rc"
        for l in gh_em gh_nvt gh_npt mh_em mh_nvt mh_npt; do
            [ -f "$dir/$l.log" ] && grep -iE 'error|fatal' "$dir/$l.log" | head -2 | sed "s/^/      $l: /"
        done
        continue
    fi
    echo "  equilibrated (EM, NVT, NPT)"
    echo "$cx READY_HOLO atoms=$atoms dry=$dry" >> "$MD/PRODUCTION_HOLO_STATUS.txt"
done

echo
echo "=== holo systems ready for production ==="
cat "$MD/PRODUCTION_HOLO_STATUS.txt" 2>/dev/null || echo "(none)"
