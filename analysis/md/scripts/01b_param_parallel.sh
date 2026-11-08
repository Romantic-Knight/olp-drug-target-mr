#!/usr/bin/env bash
# 01b_param_parallel.sh — ligand parameters, all four at once.
#
# 01_param_all.sh does the right chemistry but runs antechamber one ligand at a time, and the
# bottleneck is `sqm` (the semi-empirical AM1-BCC charge calculation), which uses barely more than
# one core. On a 32-core box that wastes most of the machine for 40-80 minutes, which is real time
# out of a 24-hour rental. Here the four ligands are parameterised concurrently instead.
#
# Everything else is identical to 01_param_all.sh: GAFF2 atom types, AM1-BCC charges, net charge
# taken from the docking manifest rather than guessed.
set -uo pipefail
cd "$(dirname "$0")/.."
MD="$(pwd)"

COMPLEXES="${COMPLEXES:-GDF11_GW8510 SUOX_Sanguinarine SUOX_Spiperone SUOX_Vorinostat}"
CHARGES="$MD/complexes/net_charges.txt"

[ -f "$CHARGES" ] || { echo "missing $CHARGES - run 01_param_all.sh once, or create it"; exit 1; }

charge_of() {
    awk -F'\t' -v w="$1" '$1==w {print $2; found=1} END{if(!found) print "NA"}' "$CHARGES"
}

echo "=== parallel ligand parameterisation (GAFF2 + AM1-BCC) ==="
pids=()
for cx in $COMPLEXES; do
    dir="$MD/complexes/$cx"
    [ -d "$dir" ] || { echo "  skip $cx (no directory)"; continue; }
    name=$(python3 -c "
import json; print(json.load(open('$MD/complexes/complex_manifest.json'))['complexes']['$cx']['ligand'])")
    q=$(charge_of "$name")
    if [ "$q" = "NA" ]; then echo "  $cx: no recorded charge for '$name'; skipping"; continue; fi

    echo "  launching $cx  ($name, charge $q)"
    (
        cd "$dir" || exit 1
        # clean any partial state from an interrupted earlier attempt
        rm -f ligand.mol2 ligand.frcmod sqm.out sqm.in ANTECHAMBER_* ATOMTYPE.INF
        antechamber -i ligand.pdb -fi pdb -o ligand.mol2 -fo mol2 \
                    -c bcc -s 2 -at gaff2 -nc "$q" -rn LIG >antechamber.log 2>&1
        rc=$?
        if [ $rc -eq 0 ] && [ -s ligand.mol2 ]; then
            parmchk2 -i ligand.mol2 -f mol2 -o ligand.frcmod -s gaff2 >parmchk2.log 2>&1
            echo "DONE $cx" > status.txt
        else
            echo "FAILED(antechamber rc=$rc) $cx" > status.txt
        fi
    ) &
    pids+=($!)
done

echo "  ${#pids[@]} ligand job(s) running; waiting..."
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=$((fail+1)); done

echo
echo "=== results ==="
ok=0
for cx in $COMPLEXES; do
    dir="$MD/complexes/$cx"
    st=$(cat "$dir/status.txt" 2>/dev/null || echo "not run")
    if [ -s "$dir/ligand.mol2" ] && [ -s "$dir/ligand.frcmod" ]; then
        atoms=$(grep -c '^ATOM' "$dir/ligand.mol2" 2>/dev/null)
        types=$(sed -n '/@<TRIPOS>ATOM/,$p' "$dir/ligand.mol2" | tail -n +2 | awk 'NF>=6{print $6}' | sort -u | tr '\n' ' ')
        printf '  OK    %-20s %3s atoms  GAFF2 types: %s\n' "$cx" "$atoms" "$types"
        ok=$((ok+1))
    else
        printf '  FAIL  %-20s %s\n' "$cx" "$st"
        tail -3 "$dir/antechamber.log" 2>/dev/null | sed 's/^/          /'
    fi
done

echo
echo "$ok/${#COMPLEXES[@]} ligands parameterised"
[ "$ok" -eq "$(echo $COMPLEXES | wc -w)" ] && echo "ALL_LIGANDS_READY"
exit 0
