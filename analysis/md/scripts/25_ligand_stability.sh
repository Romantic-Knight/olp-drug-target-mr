#!/usr/bin/env bash
# 25_ligand_stability.sh — the go/no-go check on each trajectory.
#
# Two different RMSDs answer two different questions, and conflating them is the usual mistake:
#
#   * Fitting on the **protein backbone** and measuring the ligand gives its **drift** relative to
#     the pocket. A ligand that leaves shows a monotonically rising curve here.
#   * Fitting on the **ligand itself** removes that drift and reports only **conformational**
#     change. A flat curve here with a rising one above means "it moved as a rigid body", i.e. it
#     left rather than rearranged.
#
# Also reports the minimum distance between ligand and protein over time: a ligand that has escaped
# into solvent sits far from every protein atom, which is unambiguous.
#
# Usage:  source /opt/md/activate.sh && bash 25_ligand_stability.sh [traj-dir ...]
set -uo pipefail
cd "$(dirname "$0")/.."
MD="$(pwd)"
OUT="$MD/analysis"
mkdir -p "$OUT"

dirs=("$@")
if [ ${#dirs[@]} -eq 0 ]; then
    for d in "$MD"/complexes/*/prod_rep* "$MD"/complexes/*/trial; do
        [ -s "$d/prod.tpr" ] || [ -s "$d/trial.tpr" ] && dirs+=("$d")
    done
fi

echo "=== ligand stability ==="
printf '%-30s %8s %8s %8s %8s\n' complex system "lig-drift" "lig-conf" "min-dist"
printf '%-30s %8s %8s %8s %8s\n' "" "" "(nm)" "(nm)" "(nm)"

for d in "${dirs[@]}"; do
    [ -d "$d" ] || continue
    tpr=""; xtc=""
    if [ -s "$d/prod.tpr" ]; then tpr="$d/prod.tpr"; xtc="$d/prod.xtc"
    elif [ -s "$d/trial.tpr" ]; then tpr="$d/trial.tpr"; xtc="$d/trial.xtc"; fi
    [ -n "$tpr" ] || continue
    [ -s "$xtc" ] || { echo "  $(basename $(dirname $d))/$(basename $d): no trajectory yet"; continue; }
    label="$(basename $(dirname $d))_$(basename $d)"

    # the ligand is its own residue in the topology; find its group name
    if [ ! -f "$d/index.ndx" ]; then
        printf 'Protein\nSystem\n' | gmx make_ndx -f "$tpr" -o "$d/index.ndx" >/dev/null 2>&1
    fi
    lig=$(gmx make_ndx -f "$tpr" -o /dev/null 2>/dev/null <<< "q" >/dev/null 2>&1; \
          gmx select -s "$tpr" -on /dev/null -select 'resname LIG' 2>/dev/null | head -1)
    # simplest reliable route: ask gmx for the ligand group by residue name
    liggrp=$(awk '/\[/{gsub(/[\[\]]/,"",$0); name=$0} /LIG/{print name; exit}' "$d/index.ndx" 2>/dev/null)
    if [ -z "$liggrp" ]; then
        # index groups are named after the molecule/residue; try to create one explicitly
        gmx make_ndx -f "$tpr" -o "$d/index.ndx" >/dev/null 2>&1 <<< "1 | 13
q" || true
        printf 'r LIG\nq\n' | gmx make_ndx -f "$tpr" -o "$d/index.ndx" >/dev/null 2>&1
        liggrp=$(tail -20 "$d/index.ndx" | awk '/\[/{gsub(/[\[\]]/,"",$0); n=$0} END{print n}')
    fi

    # drift: fit on protein backbone, measure the ligand
    printf 'Backbone\n'"$liggrp"'\n' 2>/dev/null | \
        gmx rms -s "$tpr" -f "$xtc" -o "$d/rmsd_lig_drift.xvg" -tu ns \
                -n "$d/index.ndx" >/dev/null 2>&1 || true
    # conformational: fit on the ligand itself
    printf "$liggrp\n$liggrp\n" | \
        gmx rms -s "$tpr" -f "$xtc" -o "$d/rmsd_lig_conf.xvg" -tu ns \
                -n "$d/index.ndx" >/dev/null 2>&1 || true
    # minimum ligand-protein distance
    printf 'Protein\n'"$liggrp"'\n' | \
        gmx mindist -s "$tpr" -f "$xtc" -od "$d/mindist.xvg" -tu ns \
                -n "$d/index.ndx" >/dev/null 2>&1 || true

    get_last() { [ -s "$1" ] && grep -v '^[#@]' "$1" | tail -1 | awk '{printf "%.3f", $2}' || echo "n/a"; }
    printf '%-30s %8s %8s %8s %8s\n' "$(basename $(dirname $d))" "$(basename $d)" \
        "$(get_last $d/rmsd_lig_drift.xvg)" "$(get_last $d/rmsd_lig_conf.xvg)" "$(get_last $d/mindist.xvg)"
done

echo
echo "READING: ligand drift < ~0.3 nm with a flat curve = the pose is stable in a pocket."
echo "         A drift that climbs past ~0.5 nm means the ligand left; do not report it as a binder."
echo "         min-dist growing past ~0.5 nm (i.e. no protein atom within 0.5 nm) means solvent exposure."
