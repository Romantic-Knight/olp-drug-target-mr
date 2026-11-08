#!/usr/bin/env bash
# 20_build_gmx.sh — Amber topology -> GROMACS system, then equilibrate. Executes, not just advises.
#
# Why Amber builds the topology and GROMACS runs the dynamics: ff19SB is an Amber force field with
# CMAP corrections and the ligand parameters came from AmberTools (GAFF2 + AM1-BCC). tleap is the
# reliable place to assemble that. `parmed` then converts to GROMACS format so production can use
# GPU GROMACS. Going the other way (pdb2gmx) would mean re-deriving the ligand parameters for a
# force field that cannot represent the metal site at all.
#
# Usage:  source /opt/md/activate.sh && bash 20_build_gmx.sh [complex ...]
set -uo pipefail
cd "$(dirname "$0")/.."
MD="$(pwd)"

COMPLEXES="${*:-GDF11_GW8510 SUOX_Sanguinarine SUOX_Spiperone SUOX_Vorinostat}"
MARGIN="${MARGIN:-1.0}"

for cx in $COMPLEXES; do
    dir="$MD/complexes/$cx"
    [ -d "$dir" ] || { echo "skip $cx"; continue; }
    [ -s "$dir/ligand.mol2" ] || { echo "skip $cx (no ligand.mol2 - parameterise first)"; continue; }
    echo "############ $cx"

    # ------------------------------------------------------------------ tleap
    # Extra ligands beyond the docked molecule. Two cases:
    #   * mte.mol2 / mte.frcmod - the SUOX molybdopterin cofactor, parameterised with GAFF2
    #     (it is a plain organic molecule; no metal). Mo itself is NOT present: AutoDock Vina
    #     cannot parse molybdenum and MCPB.py has not been run, so the metal is absent and this
    #     MUST be stated as a limitation. Keeping MTE in means the pocket geometry stays close to
    #     the one the docking pose was generated in, rather than leaving a void where the cofactor
    #     should be.
    #   * MCPB* - a bonded metal model, if MCPB.py has been run.
    metal=$(ls "$dir"/MCPB*.frcmod "$dir"/MCPB*.mol2 "$dir"/mte.frcmod "$dir"/mte.mol2 2>/dev/null \
            | xargs -r -n1 basename | tr '\n' ' ')
    if [ -n "$metal" ]; then echo "  extra parameters: $metal"; fi
    if [ -s "$dir/mte.mol2" ]; then
        echo "  NOTE: MTE included, Mo absent - record this in the paper"
    fi
    {
        echo "source leaprc.protein.ff19SB"
        echo "source leaprc.gaff2"
        echo "source leaprc.water.tip3p"
        echo "loadamberparams ligand.frcmod"
        for f in $metal; do
            case "$f" in
                *.frcmod) echo "loadamberparams $f" ;;
                *.mol2)   echo "loadmol2 $f" ;;
            esac
        done
        if [ -s "$dir/mte.mol2" ]; then echo "MTE = loadmol2 mte.mol2"; fi
        echo "LIG = loadmol2 ligand.mol2"
        echo "rec = loadpdb receptor.pdb"
        extra=""
        [ -s "$dir/mte.mol2" ] && extra="MTE"
        echo "cx  = combine { rec $extra LIG }"
        echo "check cx"
        echo "saveamberparm cx complex_dry.prmtop complex_dry.inpcrd"
        echo "solvateoct cx TIP3PBOX $MARGIN"
        echo "charge cx"
        echo "addions cx Na+ 0"
        echo "addions cx Cl- 0"
        echo "charge cx"
        echo "saveamberparm cx complex_solv.prmtop complex_solv.inpcrd"
        echo "savepdb cx complex_solv.pdb"
        echo "quit"
    } > "$dir/leap.in"

    ( cd "$dir" && tleap -f leap.in >tleap.log 2>&1 )
    if [ ! -s "$dir/complex_solv.prmtop" ]; then
        echo "  tleap FAILED"
        grep -iE 'error|fatal|could not' "$dir/tleap.log" | head -6 | sed 's/^/      /'
        continue
    fi
    atoms=$(grep -A1 '%FLAG POINTERS' "$dir/complex_solv.prmtop" | tail -1 | awk '{print $1}')
    echo "  tleap OK - solvated system: ${atoms:-?} atoms"

    # ------------------------------------------------------------------ to GROMACS
    cat > "$dir/to_gmx.py" <<'PY'
from parmed.amber import AmberParm
from parmed.gromacs import GromacsTopologyFile
p = AmberParm("complex_solv.prmtop", "complex_solv.inpcrd")
g = GromacsTopologyFile.from_structure(p)
g.save("system.top", overwrite=True)
p.save("system.gro", overwrite=True)
print("converted:", len(p.atoms), "atoms")
PY
    ( cd "$dir" && python to_gmx.py > to_gmx.log 2>&1 )
    if [ ! -s "$dir/system.top" ] || [ ! -s "$dir/system.gro" ]; then
        echo "  parmed conversion FAILED"; tail -4 "$dir/to_gmx.log" | sed 's/^/      /'; continue
    fi
    echo "  GROMACS files written ($(grep -c . "$dir/system.gro") lines)"

    # ------------------------------------------------------------------ mdp set
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
nsteps          = 50000            ; 100 ps
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
nsteps          = 50000            ; 100 ps
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

    # ------------------------------------------------------------------ EM -> NVT -> NPT
    ( cd "$dir"
      gmx grompp -f em.mdp -c system.gro -p system.top -o em.tpr -maxwarn 2 >grompp_em.log 2>&1 || exit 11
      gmx mdrun -deffnm em -nb gpu >mdrun_em.log 2>&1 || gmx mdrun -deffnm em >mdrun_em.log 2>&1 || exit 12
      gmx grompp -f nvt.mdp -c em.gro -r em.gro -p system.top -o nvt.tpr -maxwarn 2 >grompp_nvt.log 2>&1 || exit 13
      gmx mdrun -deffnm nvt -nb gpu -pme gpu -bonded gpu -update gpu >mdrun_nvt.log 2>&1 || \
      gmx mdrun -deffnm nvt >mdrun_nvt.log 2>&1 || exit 14
      gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p system.top -o npt.tpr -maxwarn 2 >grompp_npt.log 2>&1 || exit 15
      gmx mdrun -deffnm npt -nb gpu -pme gpu -bonded gpu -update gpu >mdrun_npt.log 2>&1 || \
      gmx mdrun -deffnm npt >mdrun_npt.log 2>&1 || exit 16
    )
    rc=$?
    if [ $rc -ne 0 ]; then
        echo "  equilibration FAILED at step $rc"
        for l in grompp_em grompp_nvt grompp_npt mdrun_em mdrun_nvt mdrun_npt; do
            [ -f "$dir/$l.log" ] && tail -3 "$dir/$l.log" | grep -iE 'error|fatal|fatal' | head -2 | sed "s/^/      $l: /"
        done
        continue
    fi
    echo "  equilibrated: EM -> NVT -> NPT complete"
    # a moving ligand is the go/no-go signal, so measure it right away
    ( cd "$dir"
      printf 'Protein\nProtein\n' > ndx_lig.txt 2>/dev/null
      gmx rms -s npt.tpr -f npt.xtc -o rmsd_npt.xvg -tu ns -nofit < /dev/null >/dev/null 2>&1 || true )
    echo "READY $cx"
    echo "$cx READY atoms=$atoms" >> "$MD/PRODUCTION_STATUS.txt"
done

echo
echo "=== summary ==="
cat "$MD/PRODUCTION_STATUS.txt" 2>/dev/null || echo "(nothing ready yet)"
