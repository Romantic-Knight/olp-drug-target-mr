#!/usr/bin/env bash
# 26_pbc_fix.sh — make the trajectories whole before analysing them.
#
# Why this is not optional
# ------------------------
# GROMACS writes .xtc with coordinates wrapped into the box. A ligand that drifts across a periodic
# boundary, or a protein loop that pokes out, therefore appears to teleport by roughly a box length.
# Superimposing on such coordinates is meaningless: the first analysis pass reported a protein
# backbone RMSD of ~1.8 nm over 5 ns (impossible for a folded domain) and a ligand drift of ~4 nm
# that was mostly wrapping artefact.
#
# `-pbc mol` makes each molecule whole; `-center` puts the protein in the middle of the box so the
# ligand-protein distance is not affected by where the box happens to sit; `-ur compact` packs
# everything into the smallest image which keeps the coordinates near the origin.
#
# Usage: bash 26_pbc_fix.sh [trajectory-dir ...]
#
# `set -u` must come *after* sourcing the environment: activate.sh references variables that are
# not set, so with -u already active the script died silently at the source line and produced no
# output at all. That is why the first harvest reported trajectories as "not PBC-corrected".
MD="/root/olp/analysis/md"
source /opt/md/activate.sh >/dev/null 2>&1
set -uo pipefail

if [ "$#" -gt 0 ]; then DIRS=("$@"); else DIRS=("$MD"/complexes/*/prod_rep* "$MD"/complexes/*/trial); fi

for d in "${DIRS[@]}"; do
  [ -d "$d" ] || continue
  x=$(ls "$d"/prod.xtc "$d"/trial.xtc 2>/dev/null | head -1); [ -n "$x" ] || continue
  [ -s "$x" ] || continue
  # the .tpr name varies by run: prod.tpr for production, npt.tpr/trial.tpr for the trial run,
  # so glob rather than guess
  tpr=$(ls "$d"/*.tpr 2>/dev/null | head -1); [ -n "$tpr" ] || continue
  name=$(basename "$(dirname "$d")")_$(basename "$d")
  out="$d/whole.xtc"

  # skip when the corrected file is already newer than its source
  if [ -s "$out" ] && [ "$out" -nt "$x" ]; then
    echo "  $name: up to date"
    continue
  fi

  printf 'Protein\nSystem\n' | gmx trjconv -s "$tpr" -f "$x" -o "$out" \
      -pbc mol -center -ur compact -quiet >"$d/pbc.log" 2>&1
  if [ -s "$out" ]; then
    echo "  $name: whole.xtc $(stat -c%s "$out") bytes"
  else
    echo "  $name: FAILED (see $d/pbc.log)"
  fi
done
