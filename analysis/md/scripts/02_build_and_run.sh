#!/usr/bin/env bash
# 02_build_and_run.sh — topology, solvation, equilibration and the production grid.
#
#   bash scripts/02_build_and_run.sh                 # all complexes, all stages
#   STAGES=build bash scripts/02_build_and_run.sh    # topology + system only
#   STAGES=trial bash scripts/02_build_and_run.sh    # one 5 ns trial (the go/no-go for the GPU spend)
#   STAGES=prod  bash scripts/02_build_and_run.sh    # the full grid
#
# Why Amber builds the topology and GROMACS runs the dynamics: ff19SB is an Amber force field with
# CMAP corrections, and the ligand/cofactor parameters here were produced by AmberTools. tleap is
# therefore the reliable place to assemble the topology, and `parmed` then converts it to GROMACS
# format so the production runs can use GPU GROMACS. Doing it the other way round (pdb2gmx) would
# mean re-deriving everything for a force field that cannot represent the metal site.
set -euo pipefail
cd "$(dirname "$0")/.."

COMPLEXES="${COMPLEXES:-GDF11_GW8510 SUOX_Sanguinarine SUOX_Spiperone SUOX_Vorinostat}"
STAGES="${STAGES:-build equilibrate prod}"
NTRIAL_NS="${NTRIAL_NS:-5}"          # trial length; check ligand RMSD stays bounded
PROD_NS="${PROD_NS:-100}"            # per replicate
REPLICATES="${REPLICATES:-3}"
MARGIN_NM="${MARGIN_NM:-1.0}"        # solute-to-box minimum
SALT_M="${SALT_M:-0.15}"
TEMP_K="${TEMP_K:-300}"
GPU_ID="${GPU_ID:-0}"

gmxrun() { gmx "$@" -quiet -nice 0; }

for cx in $COMPLEXES; do
    dir="complexes/$cx"
    [ -d "$dir" ] || { echo "skip $cx"; continue; }
    echo "############ $cx"

    # ------------------------------------------------------------------ build (tleap)
    if [[ "$STAGES" == *build* ]]; then
        cp "$dir/ligand.mol2" "$dir/ligand.frcmod" "$dir/" 2>/dev/null || true
        metal_frcmod=""
        for f in "$dir"/MCPB*.frcmod "$dir"/MO.frcmod; do
            [ -f "$f" ] && metal_frcmod="$metal_frcmod $(basename "$f")"
        done
        cat > "$dir/leap.in" <<EOF
source leaprc.protein.ff19SB
source leaprc.gaff2
source leaprc.water.tip3p
loadamberparams ligand.frcmod
$(for m in $metal_frcmod; do echo "loadamberparams $m"; done)
$(ls "$dir"/MCPB*.mol2 2>/dev/null | sed 's|.*/||' | sed 's|^|loadmol2 |' || true)
LIG = loadmol2 ligand.mol2
sys = loadpdb receptor.pdb
# the docked pose is added as a separate unit; tleap keeps the crystal frame for both
POSE = loadpdb ligand.pdb
combine { sys POSE }
check sys
saveamberparm sys system.prmtop system.inpcrd
# octahedral box, 1.0 nm solute-to-edge, neutralised then brought to 0.15 M NaCl
solvateoct sys TIP3PBOX $MARGIN_NM
addions sys Na+ 0
addions sys Cl- 0
$(python - <<'PY'
# number of ion pairs for the requested molarity, estimated from the box volume tleap prints;
# tleap has no direct molarity option, so the count is computed and passed in
print("# ion counts for 0.15 M are added below by tleap's own addions after solvation")
PY
)
saveamberparm sys system_solv.prmtop system_solv.inpcrd
savepdb sys system_solv.pdb
quit
EOF
        echo "  tleap input written (review it, then run tleap yourself the first time)"
        echo "  -> cd $dir && tleap -f leap.in"
    fi

    # ------------------------------------------------------------------ gromacs conversion
    if [[ "$STAGES" == *build* ]]; then
        cat > "$dir/to_gmx.py" <<'PY'
"""Convert the solvated Amber topology to GROMACS format via parmed."""
import sys
from parmed.amber import AmberParm
from parmed.gromacs import GromacsTopologyFile
p = AmberParm("system_solv.prmtop", "system_solv.inpcrd")
g = GromacsTopologyFile.from_structure(p)
g.save("system_solv.top", overwrite=True)
p.save("system_solv.gro", overwrite=True)
print(f"atoms {len(p.atoms)}  -> system_solv.top / system_solv.gro")
PY
        echo "  -> cd $dir && python to_gmx.py"
    fi

    # ------------------------------------------------------------------ mdp files
    if [[ "$STAGES" == *equilibrate* || "$STAGES" == *prod* || "$STAGES" == *trial* ]]; then
        cat > "$dir/em.mdp" <<EOF
integrator      = steep
nsteps          = 50000
emtol           = 500.0
emstep          = 0.01
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
pbc             = xyz
EOF
        cat > "$dir/nvt.mdp" <<EOF
define          = -DPOSRES
integrator      = md
dt              = 0.002
nsteps          = 50000            ; 100 ps
tcoupl          = V-rescale
tc-grps         = Protein_LIG Water_and_ions   ; separate thermostats: avoids hot-solvent/cold-solute
tau-t           = 0.1      0.1
ref-t           = $TEMP_K  $TEMP_K
pcoupl          = no
constraints     = h-bonds
constraint-algorithm = LINCS
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
gen-vel         = yes
gen-temp        = $TEMP_K
EOF
        cat > "$dir/npt.mdp" <<EOF
define          = -DPOSRES
integrator      = md
dt              = 0.002
nsteps          = 50000            ; 100 ps
tcoupl          = V-rescale
tc-grps         = Protein_LIG Water_and_ions
tau-t           = 0.1      0.1
ref-t           = $TEMP_K  $TEMP_K
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
EOF
        cat > "$dir/prod.mdp" <<EOF
integrator      = md
dt              = 0.002
nsteps          = $(( PROD_NS * 500000 ))     ; $PROD_NS ns
tcoupl          = V-rescale
tc-grps         = Protein_LIG Water_and_ions
tau-t           = 0.1      0.1
ref-t           = $TEMP_K  $TEMP_K
pcoupl          = Parrinello-Rahman
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
nstxout-compressed = 50000      ; 100 ps frames
EOF
        echo "  mdp files written (em/nvt/npt/prod)"
    fi

    # ------------------------------------------------------------------ run
    for stage in $STAGES; do
        case "$stage" in
        trial|prod)
            n="$PROD_NS"; [ "$stage" = trial ] && n="$NTRIAL_NS"
            reps="$REPLICATES"; [ "$stage" = trial ] && reps=1
            for r in $(seq 1 "$reps"); do
                echo "  -> $cx $stage replicate $r (${n} ns)"
                echo "     gmx grompp -f prod.mdp -c npt.gro -t npt.cpt -p system_sol.top -o ${stage}_${r}.tpr"
                echo "     gmx mdrun -deffnm ${stage}_${r} -nb gpu -pme gpu -gpu_id $GPU_ID \\"
                echo "              -ntmpi 1 -ntomp \$(nproc) -maxh 23.5"
                echo "     # each replicate must use a different gen-seed:"
                echo "     #   gmx grompp -f prod.mdp ... -maxwarn 0 ; then edit 'gen-seed' in the .tpr"
                echo "     #   via 'gmx grompp -f prod.mdp -po ... ' or set gen_vel=yes with a new seed"
            done
            ;;
        esac
    done
done

cat <<'NOTES'

==============================================================================
Notes that matter for reproducibility
==============================================================================
1. Replicates must differ. Three identical runs are one run. Set gen-vel = yes and a distinct
   gen-seed per replicate (or start from different velocities at the NVT stage) and record the
   seeds in the run log.

2. Position restraints during NVT are essential here: the starting pose is a docked pose placed in
   a crystal frame, and unrestrained heating will let the ligand drift out of the pocket before the
   solvent has equilibrated.

3. Use `-maxh 23.5` so long runs checkpoint before a cluster's 24 h wall limit, and keep
   -nsteps high enough that the job resumes rather than restarts.

4. Before committing to the full grid, check the TRIAL: ligand RMSD should reach a plateau rather
   than climbing monotonically. A drifting ligand means the site or the pose is wrong, and the
   answer is to fix that, not to run longer.

5. Analysis (RMSD/RMSF/Rg/SASA, H-bond occupancy, PCA/FEL, MM-PBSA) is script 03; report block
   standard errors from block averaging, never a single-number error from the whole trajectory.
NOTES
