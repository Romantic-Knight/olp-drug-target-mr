#!/usr/bin/env bash
# 21_fix_scf.sh — recover a ligand whose AM1-BCC SCF does not converge.
#
# The failure
# -----------
#     QMMM: No convergence in SCF after 1000 steps.
#     antechamber: Fatal Error! Cannot properly run sqm
#
# GW-8510 hits this: a sulfonamide sulfur plus a fused thiazole/imine/phenol system. The cause is
# the tolerance antechamber asks for, not the molecule being pathological - it writes
#
#     scfconv=1.d-10, ndiis_attempts=700,   qmcharge=0,
#
# and 1e-10 is a very tight target. Relaxing the SCF tolerance and giving DIIS more room is the
# standard, still-publishable remedy: the charges remain genuine AM1-BCC charges, only the
# convergence bookkeeping changes, and the paper should note that the SCF was extended for this
# ligand.
#
# Two earlier attempts at this were wrong and are recorded here so they are not repeated:
#   * inserting new keys into the namelist broke its structure ("Error in finding end of qmmm
#     namelist"). The existing values must be edited in place instead.
#   * treating the absence of "No convergence" as success was a false positive, because sqm had
#     died on the malformed namelist without ever running. Success must be asserted positively, by
#     finding the converged charges that antechamber needs.
#
# Usage:  source /opt/md/activate.sh && bash 21_fix_scf.sh <complex> [charge]
set -uo pipefail
cd "$(dirname "$0")/.."
MD="$(pwd)"

cx="${1:?usage: 21_fix_scf.sh <complex> [charge]}"
q="${2:-0}"
dir="$MD/complexes/$cx"
[ -d "$dir" ] || { echo "no such complex: $cx"; exit 1; }
cd "$dir"

echo "=== $cx  (net charge $q) ==="

# 1. first pass so antechamber writes sqm.in (it will fail at the sqm stage, which is expected)
if [ ! -s sqm_orig.in ]; then
    rm -f ligand.mol2 ligand.frcmod sqm.out sqm.in ANTECHAMBER_* ATOMTYPE.INF
    antechamber -i ligand_fullH.sdf -fi sdf -o ligand.mol2 -fo mol2 -c bcc -s 2 -at gaff2 \
                -nc "$q" -rn LIG >antechamber.log 2>&1 || true
fi
[ -s sqm.in ] || { echo "  no sqm.in produced; this is not an SCF failure"; tail -6 antechamber.log; exit 1; }
grep -iE 'scfconv|ndiis|maxcyc|qmcharge|qm_theory' sqm.in | sed 's/^/    original: /'

attempt=0
converged=0
for scfconv in 1.d-7 1.d-6 1.d-5; do
    attempt=$((attempt + 1))
    ndiis=$((2000 * attempt))
    # edit the values in place, leaving the namelist intact
    sed -E "s/scfconv *= *[0-9.]+d?-?[0-9]+/scfconv=${scfconv}/I; s/ndiis_attempts *= *[0-9]+/ndiis_attempts=${ndiis}/I" \
        sqm_orig.in > sqm.in
    echo "  --- attempt $attempt: scfconv=$scfconv ndiis_attempts=$ndiis"
    sqm -O -i sqm.in -o sqm.out
    # positive assertion: the converged energy AND the charges antechamber reads must be present
    if grep -qi 'No convergence' sqm.out; then
        echo "      still not converged"
        continue
    fi
    if grep -qi 'Total SCF energy' sqm.out && grep -q 'QMMM: Results' sqm.out; then
        echo "      SCF converged ($(grep -i 'Total SCF energy' sqm.out | tail -1 | sed 's/^ *//'))"
        converged=1
        break
    fi
    echo "      sqm finished without a recognisable converged result; retrying tighter"
done

if [ "$converged" -ne 1 ]; then
    echo "  SCF will not converge even at 1e-5. Do not silently substitute Gasteiger charges -"
    echo "  that would compare different charge models across complexes. Options:"
    echo "    a) a different semi-empirical Hamiltonian (qm_theory='PM3' or 'PM6')"
    echo "    b) RESP charges from a short DFT optimisation"
    echo "    c) drop this ligand from the MD set and state that in the paper"
    exit 1
fi

# 2. hand the converged charges back to antechamber
antechamber -i ligand_fullH.sdf -fi sdf -o ligand.mol2 -fo mol2 -c rc -cf sqm.out -at gaff2 \
            -nc "$q" -rn LIG >antechamber_rc.log 2>&1
if [ -s ligand.mol2 ]; then
    parmchk2 -i ligand.mol2 -f mol2 -o ligand.frcmod -s gaff2 >parmchk2.log 2>&1
    n=$(sed -n '/@<TRIPOS>ATOM/,$p' ligand.mol2 | tail -n +2 | grep -c .)
    chg=$(sed -n '/@<TRIPOS>ATOM/,$p' ligand.mol2 | tail -n +2 | awk 'NF>=9{s+=$9} END{printf "%.4f", s}')
    echo "  OK  ligand.mol2: $n atoms, summed charge $chg (requested $q)"
    echo "  frcmod revision flags: $(grep -c 'ATTN, need revision' ligand.frcmod 2>/dev/null || echo 0)"
    echo "  NOTE for the paper: AM1-BCC SCF tolerance relaxed to $(grep -oE 'scfconv=[0-9.eE+-]+' sqm.in | head -1) for this ligand"
else
    echo "  antechamber -c rc FAILED"; tail -6 antechamber_rc.log; exit 1
fi
