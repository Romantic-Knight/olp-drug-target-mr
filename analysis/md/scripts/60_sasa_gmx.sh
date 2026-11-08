#!/usr/bin/env bash
# 60_sasa_gmx.sh — solvent-accessible surface area from `gmx sasa`.
#
# Two passes, not one: `-surface`/`-output` take a single selection each, and the merged
# `Protein_LIG` group in solute.ndx is numbered against the full system, so referencing it from the
# 5725-atom solute trajectory fails ("largest allowed atom index is 5725").
# The protein-only surface is what the panel needs: the drop relative to a ligand-free reference is
# the buried interface.
#
# Deliberately no `set -e`: a single failing selection must not abort the whole sweep.
# no `set` flags at all: `set -u` before `source activate.sh` aborts the script silently,
# because that file references unset variables. Failures are handled explicitly below.
MD=/root/olp/analysis/md
OUTBASE=$MD/sasa
mkdir -p "$OUTBASE"

source /opt/md/activate.sh >/dev/null 2>&1
cd "$MD" || exit 1

echo "start $(date -Is)"
for d in complexes/*/*/; do
    d=${d%/}
    if [ ! -s "$d/solute.xtc" ] || [ ! -s "$d/solute_top.gro" ]; then
        continue
    fi
    name=$(basename "$(dirname "$d")")_$(basename "$d")

    out="$OUTBASE/${name}_protein.xvg"
    if [ ! -s "$out" ]; then
        printf 'Protein\n' | gmx sasa -s "$d/solute_top.gro" -f "$d/solute.xtc" \
            -surface Protein -output Protein -o "$out" -dt 100 > "${out%.xvg}.log" 2>&1
        if [ -s "$out" ]; then
            echo "  $name protein OK  $(grep -vc '^[#@]' "$out") rows"
        else
            echo "  $name protein FAILED"
        fi
    else
        echo "  $name protein cached"
    fi
done
echo "end $(date -Is)"
