#!/usr/bin/env bash
# 70_fpocket.sh — pocket detection on the two receptors, for the druggability panel.
#
# What this produces, per receptor: the number of alpha spheres, the pocket volume and the
# druggability score of the top-ranked pocket. fpocket ranks cavities by a score that already
# weights volume, hydrophobicity and enclosure, so the top pocket is the right one to compare.
#
# Both receptors come from the MD equilibration, stripped of solvent and cofactor, so the
# comparison is like-for-like: same force field, same treatment.
set -o pipefail
MD=/root/olp/analysis/md
OUT=$MD/fpocket
mkdir -p "$OUT"

one () {  # $1 = label, $2 = pdb, $3 = which chain/selection description
    local label="$1" pdb="$2"
    local d="$OUT/$label"
    rm -rf "$d"; mkdir -p "$d"
    cp "$pdb" "$d/receptor.pdb"
    ( cd "$d" && fpocket -f receptor.pdb > fpocket.log 2>&1 )
    local info
    info=$(ls "$d"/receptor_out/*_info.txt 2>/dev/null | head -1)
    if [ -z "$info" ]; then
        echo "  $label: fpocket produced no info file"
        tail -4 "$d/fpocket.log" | sed 's/^/      /'
        return 1
    fi
    # the info file lists Pocket 1..N with druggability score, volume, etc.
    echo "  $label:"
    awk '/Pocket [0-9]+ :/ {p=$2} /Druggability Score/ {ds=$4} /Volume/ {v=$2}
         /Pocket [0-9]+ :/ {n++}
         n<=3 && /Druggability Score|Volume|Number of Alpha Spheres/ {print "      pocket " p ": " $0}
         n==1 && /^Pocket 1/ {print "      --- top pocket ---"}
        ' "$info" | head -12
    # also write a tidy summary line
    python3 - "$info" "$label" <<'PY'
import re, sys
info, label = sys.argv[1], sys.argv[2]
txt = open(info).read()
blocks = re.split(r"\nPocket (\d+) :\n", txt)
rows = []
for i in range(1, len(blocks), 2):
    pid = int(blocks[i]); body = blocks[i+1]
    def g(pat, cast=float):
        m = re.search(pat + r"\s*:\s*([-\d.eE+]+)", body)
        return cast(m.group(1)) if m else None
    rows.append(dict(pocket=pid,
                     druggability=g(r"Druggability Score"),
                     volume=g(r"Volume"),
                     alpha_spheres=g(r"Number of Alpha Spheres", int),
                     mean_alpha=g(r"Mean alpha sphere radius"),
                     hydrophobicity=g(r"Hydrophobicity score"),
                     polar=g(r"Polarity score")))
rows = [r for r in rows if r["druggability"] is not None]
rows.sort(key=lambda r: -r["druggability"])
import json
json.dump(dict(label=label, n_pockets=len(rows), top=rows[0] if rows else None,
               top5=rows[:5]),
          open(f"/root/olp/analysis/md/fpocket/{label}.json", "w"), indent=2)
print(f"      -> {len(rows)} pockets, top druggability {rows[0]['druggability']:.3f}" if rows else "      -> no scored pockets")
PY
}

echo "start $(date -Is)"
one SUOX   "$MD/fpocket/inputs/SUOX.pdb"
one GDF11  "$MD/fpocket/inputs/GDF11.pdb"
echo "end $(date -Is)"
