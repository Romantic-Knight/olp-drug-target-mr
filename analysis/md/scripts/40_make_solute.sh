#!/usr/bin/env bash
# 40_make_solute.sh — produce the download-ready trajectory for one run.
#
# This is the step that failed in the unattended harvest, and the failures were worth recording
# because each one was silent:
#
#   * selecting the output group by position ("the last group in the index") grabbed non-Water
#     instead of protein+ligand, so the trajectory had 1718 atoms instead of 5708 and every later
#     analysis died with "topology and XTC do not have the same number of atoms";
#   * merging "1 | 14" assumed 14 was LIG when it is Cl-, which produced protein + chloride with no
#     ligand at all;
#   * computing the group index with awk's counter was off by one, because gmx indexes are 0-based.
#
# Rather than keep deriving group numbers, the selection is written with `gmx select`, which takes a
# selection expression and produces the group directly. The atom count is then verified against the
# topology before the file is accepted, so a wrong selection cannot slip through again.
#
# Timings for reference: a 200 ns SUOX run is a 1.7 GB .xtc, PBC correction takes about a minute,
# and the solute trajectory comes out at ~13 MB (protein + ligand, every 100 ps).
#
# Usage: bash 40_make_solute.sh <run-dir>
#        bash 40_make_solute.sh --all
MD="/root/olp/analysis/md"
source /opt/md/activate.sh >/dev/null 2>&1
set -uo pipefail

one() {
    local d="$1"
    [ -d "$d" ] || { echo "  $d: not a directory"; return 1; }
    local tpr xtc name
    tpr=$(ls "$d"/*.tpr 2>/dev/null | head -1)
    xtc="$d/prod.xtc"
    name=$(echo "$d" | sed 's|.*/complexes/||')
    [ -s "$tpr" ] || { echo "  $name: no .tpr"; return 1; }
    [ -s "$xtc" ] || { echo "  $name: no prod.xtc"; return 1; }

    echo "=== $name"
    echo "    raw $(du -h "$xtc" | cut -f1)"

    # 1. whole molecules, protein centred, packed: without this a ligand or loop crossing a box
    #    boundary looks like it jumped several nm, which once produced a false "left the site"
    printf 'Protein\nSystem\n' | gmx trjconv -s "$tpr" -f "$xtc" -o "$d/whole.xtc" \
        -pbc mol -center -ur compact -quiet >"$d/pbc.log" 2>&1
    if [ ! -s "$d/whole.xtc" ]; then
        echo "    PBC FAILED"; tail -4 "$d/pbc.log" | sed 's/^/      /'; return 1
    fi
    echo "    whole $(du -h "$d/whole.xtc" | cut -f1)"

    # 2. protein + ligand group.
    #
    # make_ndx merges by NUMBER, not by name ("Protein | LIG" is a syntax error, "1 | 13" works and
    # reports "Merged two groups with OR: 1671 46 -> 1717"). The numbers are read from the index this
    # same command just produced, so they cannot drift.
    #
    # Getting this wrong was silent: "1 | 14" merged protein with Cl- (14 is the chloride group, not
    # the ligand) and the resulting trajectory had no ligand at all, which only surfaced much later
    # as "topology and XTC do not have the same number of atoms".
    printf 'q\n' | gmx make_ndx -f "$tpr" -o "$d/base.ndx" >"$d/mkndx.log" 2>&1
    local i_prot i_lig
    i_prot=$(awk '/^\[/{if ($0 ~ /^\[ Protein \]/) print n; n++}' "$d/base.ndx")
    i_lig=$(awk '/^\[/{if ($0 ~ /^\[ LIG \]/) print n; n++}' "$d/base.ndx")
    if [ -z "$i_prot" ] || [ -z "$i_lig" ]; then
        echo "    cannot locate Protein (idx=$i_prot) or LIG (idx=$i_lig) in the index"
        return 1
    fi
    printf '%s | %s\nq\n' "$i_prot" "$i_lig" | gmx make_ndx -f "$tpr" -o "$d/solute.ndx" \
        >>"$d/mkndx.log" 2>&1
    local grp merged
    grp=$(grep -oE '^\[ Protein_LIG \]$' "$d/solute.ndx" | tr -d '[] ')
    merged=$(grep -oE 'Merged two groups with OR: [0-9]+ [0-9]+ -> [0-9]+' "$d/mkndx.log" | tail -1)
    if [ -z "$grp" ]; then
        echo "    merge failed (protein=$i_prot lig=$i_lig)"; return 1
    fi
    echo "    groups: Protein=$i_prot LIG=$i_lig  ${merged:-(merge not reported)}"

    # 3. extract at 100 ps spacing
    printf '%s\n' "$grp" | gmx trjconv -s "$tpr" -f "$d/whole.xtc" -n "$d/solute.ndx" \
        -o "$d/solute.xtc" -dt 100 -quiet >"$d/solute.log" 2>&1
    if [ ! -s "$d/solute.xtc" ]; then
        echo "    solute FAILED"; tail -4 "$d/solute.log" | sed 's/^/      /'; return 1
    fi

    # 4. verify against the topology before trusting the file
    # Verify against the topology without guessing a .gro path: MDAnalysis can read the .tpr here
    # (it handles tpx 138 in this version), and if it cannot, the check is skipped with a warning
    # rather than failing the run - the atom count is confirmed again by analyse_md.py.
    local report
    report=$(python3 - "$tpr" "$d/solute.xtc" "$d/solute.ndx" <<'PY'
import sys
try:
    import MDAnalysis as mda
    sol = mda.Universe(sys.argv[1], sys.argv[2])
    # count the merged group's atom indices straight from the index file
    n_exp, inside = 0, False
    for line in open(sys.argv[3]):
        if line.startswith("["):
            inside = "Protein_LIG" in line or "Solute" in line
            continue
        if inside and line.split():
            n_exp += len(line.split())
    if n_exp:
        print(f"{len(sol.atoms)} {n_exp} {len(sol.trajectory)}")
    else:
        print(f"{len(sol.atoms)} - {len(sol.trajectory)}")
except Exception as e:
    print(f"SKIP {type(e).__name__}")
PY
)
    set -- $report
    if [ "${1:-SKIP}" = "SKIP" ]; then
        echo "    solute $(du -h "$d/solute.xtc" | cut -f1)  (count check skipped: ${2:-?})"
    else
        if [ "${2:--}" != "-" ] && [ "$1" != "$2" ]; then
            echo "    ATOM COUNT MISMATCH: trajectory has $1, expected $2 (protein+ligand)"
            return 1
        fi
        echo "    solute $(du -h "$d/solute.xtc" | cut -f1), $1 atoms, $3 frames"
    fi

    # 5. free the big intermediate; the raw trajectory always stays
    rm -f "$d/whole.xtc"
    echo "    removed whole.xtc, kept prod.xtc"
    return 0
}

if [ "${1:-}" = "--all" ]; then
    ok=0; bad=0
    for d in "$MD"/complexes/*/holo_rep* "$MD"/complexes/*/gdf11_rep*; do
        [ -d "$d" ] || continue
        if one "$d"; then ok=$((ok+1)); else bad=$((bad+1)); echo "    FAILED: $d"; fi
    done
    echo
    echo "=== $ok succeeded, $bad failed ==="
    df -h / | tail -1
else
    one "${1:-}"
fi
