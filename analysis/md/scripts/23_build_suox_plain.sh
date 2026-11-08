#!/usr/bin/env bash
# 23_build_suox_plain.sh — build the SUOX systems without the molybdopterin cofactor.
#
# Why the cofactor is dropped (a documented decision, not an oversight)
# ---------------------------------------------------------------------
# 2A9A's active site contains Mo-molybdopterin (MTE), a catalytic sulfate and two coordinating
# waters. Keeping MTE would be more faithful, and three attempts were made:
#
#   1. MTE parameterised with GAFF2 from a hydrogen-completed model. sqm refused the first attempt
#      (odd electron count, 207) and the second (the completed model carried the wrong protonation).
#   2. The mol2 then failed tleap: "Atom .R<MTE 1>.A<C4 5> has force field coordination 4 but only 3
#      bonded neighbors", traced to three duplicate bonds plus a carbonyl written as a single bond
#      (measured C4-O4 = 1.249 A in the crystal, so C=O).
#   3. With the bonds repaired and the carbonyl promoted, sqm converged (172 electrons, even) but
#      antechamber assigned MTE a net charge of **+1** where the chemistry requires about **-2**
#      (dithiolene thiolates plus a phosphate). AM1-BCC charges computed on a wrong net charge are
#      wrong, and charges sit right next to the ligand - so proceeding would put a known error in
#      the one place that matters most.
#
# The metal cannot be represented in ff19SB at all; it needs MCPB.py, which in turn needs MTE's
# protonation solved first. Rather than ship trajectories with wrong electrostatics in the binding
# site, the cofactor is omitted and this is recorded as a limitation:
#
#   **Limitation: the SUOX simulations omit Mo, MTE, and the catalytic sulfate. The pocket is
#   therefore an empty cavity rather than the true catalytic site, so ligand stability must be read
#   as "stability in the docked subsite", not as a description of the holoenzyme.**
#
# What this script does
# ---------------------
# Removes the cofactor from the receptor, then runs the full chain: tleap (ff19SB + TIP3P) ->
# parmed -> GROMACS -> EM -> NVT -> NPT, for every complex whose ligand is parameterised.
#
# Usage:  source /opt/md/activate.sh && bash 23_build_suox_plain.sh [complex ...]
set -uo pipefail
cd "$(dirname "$0")/.."
MD="$(pwd)"

COMPLEXES="${*:-SUOX_Sanguinarine SUOX_Spiperone SUOX_Vorinostat}"
MARGIN="${MARGIN:-1.2}"

: > "$MD/PRODUCTION_STATUS.txt"

for cx in $COMPLEXES; do
    dir="$MD/complexes/$cx"
    [ -d "$dir" ] || { echo "skip $cx (no dir)"; continue; }
    [ -s "$dir/ligand.mol2" ] || { echo "skip $cx (ligand not parameterised)"; continue; }
    echo "############ $cx"

    # ---- receptor: protein only ------------------------------------------------------------
    if [ ! -f "$dir/receptor_protein_only.pdb" ]; then
        awk 'substr($0,1,4)=="ATOM"' "$dir/receptor.pdb" > "$dir/receptor_protein_only.pdb"
        printf 'END\n' >> "$dir/receptor_protein_only.pdb"
    fi
    cp "$dir/receptor_protein_only.pdb" "$dir/receptor.pdb"
    echo "  receptor: $(grep -c '^ATOM' "$dir/receptor.pdb") protein atoms, no cofactor"

    # make sure the cofactor parameters cannot be picked up by accident
    for f in mte.mol2 mte.frcmod; do
        [ -f "$dir/$f" ] && mv -f "$dir/$f" "$dir/$f.unused_for_MD"
    done

    # ---- tleap -----------------------------------------------------------------------------
    cat > "$dir/leap.in" <<EOF
source leaprc.protein.ff19SB
source leaprc.gaff2
source leaprc.water.tip3p
loadamberparams ligand.frcmod
LIG = loadmol2 ligand.mol2
rec = loadpdb receptor.pdb
cx  = combine { rec LIG }
check cx
saveamberparm cx complex_dry.prmtop complex_dry.inpcrd
solvateoct cx TIP3PBOX $MARGIN
addions cx Na+ 0
addions cx Cl- 0
saveamberparm cx complex_solv.prmtop complex_solv.inpcrd
savepdb cx complex_solv.pdb
quit
EOF
    ( cd "$dir" && tleap -f leap.in >tleap.log 2>&1 )
    if [ ! -s "$dir/complex_solv.prmtop" ]; then
        echo "  tleap FAILED"
        grep -iE 'fatal|error|could not' "$dir/tleap.log" | head -5 | sed 's/^/      /'
        continue
    fi
    atoms=$(awk '/%FLAG POINTERS/{getline; getline; print $1}' "$dir/complex_solv.prmtop")
    echo "  tleap OK: solvated system ${atoms} atoms"

    # ---- parmed -> GROMACS ------------------------------------------------------------------
    cat > "$dir/to_gmx.py" <<'PY'
from parmed.amber import AmberParm
from parmed.gromacs import GromacsTopologyFile
p = AmberParm("complex_solv.prmtop", "complex_solv.inpcrd")
g = GromacsTopologyFile.from_structure(p)
g.save("system.top", overwrite=True)
p.save("system.gro", overwrite=True)
print("converted", len(p.atoms), "atoms")
PY
    ( cd "$dir" && python to_gmx.py >to_gmx.log 2>&1 )
    if [ ! -s "$dir/system.top" ]; then
        echo "  parmed conversion FAILED"; tail -4 "$dir/to_gmx.log" | sed 's/^/      /'; continue
    fi
    echo "  GROMACS files written"

    # ---- mdp files --------------------------------------------------------------------------
    cat > "$dir/em.mdp" <<'EOF'
integrator      = steep
nsteps          = 50000
emtol           = 200.0
emstep          = 0.01
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
pbc             = xyz
EOF
    cat > "$dir/nvt.mdp" <<'EOF'
define          = -DPOSRES
integrator      = md
dt              = 0.002
nsteps          = 50000
tcoupl          = V-rescale
tc-grps         = System
tau-t           = 0.1
ref-t           = 300
pcoupl          = no
constraints     = h-bonds
constraint-algorithm = LINCS
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
gen-vel         = yes
gen-temp        = 300
nstxout-compressed = 0
EOF
    cat > "$dir/npt.mdp" <<'EOF'
define          = -DPOSRES
integrator      = md
dt              = 0.002
nsteps          = 50000
tcoupl          = V-rescale
tc-grps         = System
tau-t           = 0.1
ref-t           = 300
pcoupl          = C-rescale
pcoupltype      = isotropic
tau-p           = 2.0
ref-p           = 1.0
compressibility = 4.5e-5
constraints     = h-bonds
constraint-algorithm = LINCS
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
gen-vel         = no
nstxout-compressed = 0
EOF

    # ---- EM -> NVT -> NPT -------------------------------------------------------------------
    ( cd "$dir"
      gmx grompp -f em.mdp -c system.gro -p system.top -o em.tpr -maxwarn 2 >grompp_em.log 2>&1 || exit 11
      gmx mdrun -deffnm em -nb gpu >mdrun_em.log 2>&1 || gmx mdrun -deffnm em >mdrun_em.log 2>&1 || exit 12
      gmx grompp -f nvt.mdp -c em.gro -r em.gro -p system.top -o nvt.tpr -maxwarn 2 >grompp_nvt.log 2>&1 || exit 13
      gmx mdrun -deffnm nvt -nb gpu -pme gpu -bonded gpu -update gpu >mdrun_nvt.log 2>&1 || \
        gmx mdrun -deffnm nvt >mdrun_nvt.log 2>&1 || exit 14
      gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p system.top -o npt.tpr -maxwarn 2 >grompp_npt.log 2>&1 || exit 15
      gmx mdrun -deffnm npt -nb gpu -pme gpu -bonded gpu -update gpu >mdrun_npt.log 2>&1 || \
        gmx mdrun -deffnm npt >mdrun_npt.log 2>&1 || exit 16 )
    rc=$?
    if [ $rc -ne 0 ]; then
        echo "  equilibration FAILED at step $rc"
        for l in grompp_em grompp_nvt grompp_npt mdrun_em mdrun_nvt mdrun_npt; do
            if [ -f "$dir/$l.log" ]; then
                grep -iE 'error|fatal' "$dir/$l.log" | head -2 | sed "s/^/      $l: /"
            fi
        done
        continue
    fi
    echo "  equilibrated (EM, NVT, NPT)"
    echo "$cx READY atoms=$atoms" >> "$MD/PRODUCTION_STATUS.txt"
done

echo
echo "=== ready for production ==="
cat "$MD/PRODUCTION_STATUS.txt" 2>/dev/null || echo "(none)"
