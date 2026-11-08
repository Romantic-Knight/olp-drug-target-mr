#!/usr/bin/env bash
# 41_analyse_all.sh — stability analysis for every run, using a matching topology.
#
# The trap this avoids: analyse_md.py pairs the trajectory with the system .gro, but that .gro is the
# solvated system (e.g. 14,579 atoms) while solute.xtc holds only protein+ligand (1,717). MDAnalysis
# then raises "topology and XTC do not have the same number of atoms" — which is exactly the error
# that made the first harvest report nothing for the holo series.
#
# Fix: dump frame 0 of solute.xtc to a .gro, and use THAT as the topology. It is by construction the
# same atom set and order as the trajectory.
MD=/root/olp/analysis/md
source /opt/md/activate.sh >/dev/null 2>&1
cd "$MD" || exit 1

for d in complexes/*/*/; do
    d=${d%/}
    [ -s "$d/solute.xtc" ] || continue
    [ -s "$d/solute.ndx" ] || continue

    # Topology that matches the solute trajectory.
    #
    # Three traps here, all of which produce "Index[N] ... is larger than the number of atoms":
    #   * solute.ndx is numbered against the FULL system (protein 1-5669, ligand 5715-5753, because
    #     MTE/SO4 occupy 5670-5714), while solute.xtc holds only the 5708 selected atoms, so the
    #     index cannot be re-applied to the subset;
    #   * selecting group 0 of the 23,779-atom prod.tpr against a 5,708-atom trajectory cannot work
    #     either, because the two atom sets differ;
    #   * so the only consistent source is the FULL trajectory (prod.xtc) with the same index and tpr.
    # Dumping frame 0 from prod.xtc through that index yields exactly the solute atom set and order.
    if [ ! -s "$d/solute_top.gro" ]; then
        grp=$(grep -oE '^\[ (Protein_LIG|Solute) \]' "$d/solute.ndx" | head -1 | tr -d '[]')
        [ -z "$grp" ] && grp=Protein_LIG
        printf '%s\n' "$grp" | gmx trjconv -s "$d/prod.tpr" -f "$d/prod.xtc" -n "$d/solute.ndx" \
            -o "$d/solute_top.gro" -dump 0 >"$d/gro_dump.log" 2>&1
        if [ -s "$d/solute_top.gro" ]; then
            echo "    topology: $(awk 'NR==2{print $1}' "$d/solute_top.gro") atoms (from $grp)"
        else
            echo "    topology FAILED:"; grep -A3 Fatal "$d/gro_dump.log" | sed 's/^/      /'; continue
        fi
    fi

done

# One invocation for ALL runs: analyse_md.py writes a single md_stability.json, so calling it per
# directory silently keeps only the last run (that is how the first harvest ended up reporting one
# record). With no arguments it globs every run itself.
echo
echo "=== stability analysis (single pass over every run) ==="
python3 scripts/analyse_md.py 2>&1 | grep -E '^[A-Za-z]|left the site|stable|intermediate' | sed 's/^/  /'
echo
echo "  -> $(python3 -c "import json;d=json.load(open('analysis/md_stability.json'));print(len(d),'records')" 2>/dev/null)"
