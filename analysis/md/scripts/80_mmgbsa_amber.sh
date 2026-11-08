#!/usr/bin/env bash
# 80_mmgbsa_amber.sh — MM-GBSA with Amber's own MMPBSA.py.
#
# Why not gmx_MMPBSA: it crashes in res2map on these systems (Amber->GROMACS converted topologies
# with hand-built indices). Amber's MMPBSA.py takes the original prmtop and a cpptraj trajectory,
# so there is no residue remapping step to go wrong.
#
# Inputs per run:
#   <complex>/complex_dry.prmtop   topology of protein + ligand, no solvent (5708 / 1717 atoms)
#   <run>/solute.xtc               the analysed trajectory (same atom count, PBC-corrected)
# Output: <run>/mmgbsa/FINAL_RESULTS_MMPBSA.dat
#
# 100 evenly spaced frames is enough for a converged mean (single-trajectory MM-GBSA) and keeps
# the runtime at a few minutes per run on 32 cores.
MD=/root/olp/analysis/md
FRAMES=${FRAMES:-100}
INTERVAL=${INTERVAL:-33}   # 2001 frames / 33 -> ~60 samples across the run

one () {
    local run="$1" cx="$2"
    local d="$MD/complexes/$cx/$run"
    local prm="$MD/complexes/$cx/mmpbsa.prmtop"
    [ -s "$prm" ] || { echo "  $run: no mmpbsa.prmtop"; return 1; }
    [ -s "$d/solute.xtc" ] || { echo "  $run: no solute.xtc"; return 1; }

    local natom_prm natom_xtc
    natom_prm=$(awk '/^%FLAG POINTERS/{getline; getline; print $1}' "$prm")
    natom_gro=$(awk 'NR==2{print $1}' "$d/solute_top.gro" 2>/dev/null)
    if [ -n "$natom_gro" ] && [ "$natom_prm" != "$natom_gro" ]; then
        echo "  $run: SKIP - prmtop $natom_prm vs trajectory $natom_gro atoms"
        return 1
    fi
    echo "  $run: prmtop $natom_prm atoms (matches trajectory)"

    mkdir -p "$d/mmgbsa"
    cd "$d/mmgbsa" || return 1

    # GROMACS -> Amber trajectory, using the prmtop as the topology so atom order is guaranteed.
    # receptor and ligand topologies: MMPBSA.py only subtracts them if it is given both, so with
    # -cp alone it reports the complex total and no DELTA. Split them from the complex prmtop.
    cat > split.in <<EOF
parm $prm
parmstrip :LIG
parmwrite out rec.prmtop
clear all
parm $prm
parmstrip !:LIG
parmwrite out lig.prmtop
go
EOF
    cpptraj -i split.in > split.log 2>&1 || true

    cat > cpptraj.in <<EOF
parm $prm
trajin $d/solute.xtc
trajout traj.mdcrd
go
EOF
    cpptraj -i cpptraj.in > cpptraj.log 2>&1
    if [ ! -s traj.mdcrd ]; then
        echo "  $run: cpptraj conversion FAILED"; tail -5 cpptraj.log | sed 's/^/      /'; return 1
    fi

    cat > mmpbsa.in <<EOF
&general
   startframe=1, endframe=100000, interval=$INTERVAL,
   keep_files=0, verbose=1,
/
&gb
   igb=5, saltcon=0.150,
/
EOF

    # a single topology is enough: MMPBSA.py derives receptor and ligand from the mask
    local extra=""
    [ -s rec.prmtop ] && [ -s lig.prmtop ] && extra="-rp rec.prmtop -lp lig.prmtop"
    MMPBSA.py -O -i mmpbsa.in -o FINAL_RESULTS_MMPBSA.dat \
        -cp "$prm" $extra -y traj.mdcrd > mmpbsa.log 2>&1

    if [ -s FINAL_RESULTS_MMPBSA.dat ]; then
        local dg
        dg=$(awk '/^DELTA TOTAL/{print $3}' FINAL_RESULTS_MMPBSA.dat | tail -1)
        echo "  $run: DELTA TOTAL = ${dg} kcal/mol"
    else
        echo "  $run: MMPBSA FAILED"; tail -6 mmpbsa.log | sed 's/^/      /'; return 1
    fi
}

echo "start $(date -Is)  frames=$FRAMES"
if [ -n "${1:-}" ]; then
    one "$1" "$2"
else
    for cx in SUOX_Sanguinarine SUOX_Spiperone SUOX_Vorinostat GDF11_GW8510; do
        for d in "$MD/complexes/$cx"/holo_rep* "$MD/complexes/$cx"/gdf11_rep*; do
            [ -d "$d" ] || continue
            one "$(basename "$d")" "$cx"
        done
    done
fi
echo "end $(date -Is)"
