#!/usr/bin/env bash
# 50_mmpbsa.sh — binding free energy (MM-PBSA) for one run.
#
# Input requirements, learned the hard way:
#   * gmx_MMPBSA refuses solvent: "-cs" must be a DRY complex (protein+ligand), not the solvated tpr
#     ("Complex included water or ions"). It also refuses Amber prmtop via -cp (that flag takes a
#     GROMACS .tpr/.top); the earlier attempts failed for exactly those two reasons.
#   * So we make a dry tpr from the dry topology, pair it with the PBC-corrected solute trajectory,
#     and hand the same index that maps group 1 = receptor, group 2 = ligand.
#
# Usage: 50_mmpbsa.sh <run-dir>   (e.g. complexes/SUOX_Spiperone/holo_rep1)
set -uo pipefail

d=${1:?usage: 50_mmpbsa.sh <run-dir>}
source /opt/md/activate.sh >/dev/null 2>&1
cd "$d" || exit 1

echo "=== $(basename "$(dirname "$d")")/$(basename "$d") ==="

# 1. dry tpr: pull protein+ligand out of the solvated tpr with the same index used for solute.xtc
if [ ! -s dry.tpr ]; then
    printf 'Protein_LIG\n' | gmx convert-tpr -s prod.tpr -n solute.ndx -o dry.tpr \
        >dry_tpr.log 2>&1
    if [ ! -s dry.tpr ]; then
        echo "  dry.tpr FAILED"; grep -A3 -i error dry_tpr.log | head -6; exit 1
    fi
fi
echo "  dry.tpr: $(stat -c%s dry.tpr) bytes"

# 2. index for gmx_MMPBSA: group 1 = receptor (protein), group 2 = ligand
if [ ! -s mmpbsa.ndx ]; then
    printf '1\n13\nq\n' | gmx make_ndx -f dry.tpr -o mmpbsa.ndx >mmpbsa_ndx.log 2>&1
    echo "  mmpbsa.ndx groups: $(grep -c '^\[' mmpbsa.ndx 2>/dev/null)"
fi

# 3. sanity: does the dry tpr match the trajectory atom count?
n_tpr=$(gmx check -s1 dry.tpr 2>&1 | grep -oE '# Atoms *[0-9]+' | grep -oE '[0-9]+' | head -1)
n_xtc=$(gmx check -f solute.xtc 2>&1 | grep -oE '# Atoms *[0-9]+' | grep -oE '[0-9]+' | head -1)
echo "  atoms: dry.tpr=$n_tpr  solute.xtc=$n_xtc"
if [ -n "$n_tpr" ] && [ -n "$n_xtc" ] && [ "$n_tpr" != "$n_xtc" ]; then
    echo "  ATOM MISMATCH — stopping"; exit 1
fi

# 4. run MM-PBSA on a strided subset first (a full 2001-frame run is hours; 100 frames is enough
#    for a mean and keeps the iteration fast while we validate the setup)
if [ ! -s FINAL_RESULTS_MMPBSA.dat ] || [ "${FORCE:-0}" = "1" ]; then
    gmx_MMPBSA -O -i mmpbsa.in -cs dry.tpr -ci mmpbsa.ndx -cg 1 2 \
        -ct solute.xtc -cp dry.tpr -o FINAL_RESULTS_MMPBSA.dat \
        -eo FINAL_RESULTS_MMPBSA.csv -nogui > gmx_mmpbsa.log 2>&1
fi

if [ -s FINAL_RESULTS_MMPBSA.dat ]; then
    grep -E 'DELTA TOTAL|DELTA G binding|GGAS|GSOLV' FINAL_RESULTS_MMPBSA.dat | head -6 | sed 's/^/    /'
else
    echo "  MM-PBSA FAILED:"; tail -12 gmx_mmpbsa.log | sed 's/^/    /'
fi
