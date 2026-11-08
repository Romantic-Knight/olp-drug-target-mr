#!/usr/bin/env bash
# 51_mmpbsa_clean.sh — MM-PBSA with a two-group index, built the way gmx_MMPBSA expects.
#
# What went wrong repeatedly:
#   * "make_ndx -f dry.tpr -o mmpbsa.ndx" COPIES every default group first (System, Protein,
#     Protein-H, ... LIG), and only then appends the ones you asked for. So -cg 1 2 pointed at
#     Protein-H, not the ligand, and gmx_MMPBSA crashed in res2map with
#     "string indices must be integers" because it was trying to treat 2863 protein atoms as one
#     ligand residue.
#   * -cs and -cp cannot both be the same file (duplicate-argument error), and -cp must be a GAFF
#     .top, not a .tpr.
#
# This script therefore builds an index that contains EXACTLY two groups, in the order
# receptor-then-ligand, and passes -cg 1 2.
set -uo pipefail

d=${1:?usage: 51_mmpbsa_clean.sh <run-dir>}
FRAMES=${FRAMES:-100}
source /opt/md/activate.sh >/dev/null 2>&1
cd "$d" || exit 1
echo "=== $(basename "$(dirname "$d")")/$(basename "$d")  (${FRAMES} frames) ==="

[ -s dry.tpr ] || { echo "  no dry.tpr"; exit 1; }
[ -s dry.top ] || { echo "  no dry.top"; exit 1; }

# 1. two-group index, written from scratch (gmx select honours the order given)
rm -f two.ndx
gmx select -s dry.tpr -select 'group "Protein"' -on two.ndx >sel1.log 2>&1
gmx select -s dry.tpr -select 'group "LIG"'     -on two.ndx -on two.ndx >sel2.log 2>&1 || true

# the reliable way: hand-write the two blocks from the default index
printf 'q\n' | gmx make_ndx -f dry.tpr -o all.ndx >/dev/null 2>&1
python3 - <<'PY'
import re
txt = open("all.ndx").read()
blocks = re.findall(r"\[ (.+?) \]\n((?:[^\[].*\n)*)", txt)
want = {}
for name, body in blocks:
    if name not in want:
        want[name] = body
out = []
for name in ("Protein", "LIG"):
    if name in want:
        out.append(f"[ {name} ]\n{want[name]}")
open("two.ndx", "w").write("".join(out))
print("   two.ndx groups:", [n for n in ("Protein", "LIG") if n in want])
PY
echo "  two.ndx: $(grep -c '^\[' two.ndx) groups, $(awk 'NF&&!/^\[/{n+=NF} END{print n}' two.ndx) atoms"

# 2. input file: strided frames only, GB implicit solvent
last=$(gmx check -f solute.xtc 2>&1 | awk '/^Coords/{print $2}')
int=$(( last / FRAMES )); [ "$int" -lt 1 ] && int=1
cat > mmpbsa.in <<EOF
&general
  startframe=1, endframe=${last}, interval=${int},
  verbose=2,
/
&pb
  istrng=0.150,
/
EOF
echo "  frames ${last}, interval ${int}"

# 3. run
rm -f FINAL_RESULTS_MMPBSA.dat
gmx_MMPBSA -O -i mmpbsa.in -ci two.ndx -cg 1 2 \
    -cs dry.tpr -ct solute.xtc -cp dry.top \
    -o FINAL_RESULTS_MMPBSA.dat -nogui > gmx_MMPBSA.log 2>&1

if [ -s FINAL_RESULTS_MMPBSA.dat ]; then
    echo "  OK"
    sed -n '/DELTA TOTAL/,+2p' FINAL_RESULTS_MMPBSA.dat | head -4 | sed 's/^/    /'
    grep -E '^DELTA G binding' FINAL_RESULTS_MMPBSA.dat | head -2 | sed 's/^/    /'
else
    echo "  FAILED — last lines:"
    tail -6 gmx_MMPBSA.log | sed 's/^/    /'
fi
