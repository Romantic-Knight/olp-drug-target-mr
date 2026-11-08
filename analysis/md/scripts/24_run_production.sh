#!/usr/bin/env bash
# 24_run_production.sh — the production grid.
#
# 3 complexes x 3 independent replicates, run CONCURRENTLY (the benchmark showed one 4090 is not
# saturated by a ~25k atom system, so a few simultaneous mdrun processes give more aggregate ns/day
# than one at a time: N=2 gave 979 ns/day aggregate vs 768 for N=1).
#
# Design points that matter scientifically
# ----------------------------------------
# * **Three replicates with distinct velocity seeds.** Three identical runs are one run; the whole
#   point of replicates is to show the binding mode is reproducible rather than a single lucky
#   trajectory. gen-seed is set per replicate and recorded.
# * **Checkpointing with -maxh** so a run can resume if the session or the box restarts, and so
#   nothing is lost when the rental window closes.
# * **.xtc output every 10 ps**, compressed. Never .trr - it would fill the disk.
# * **No position restraints**, unlike the equilibration. Restraints exist to stop the ligand being
#   thrown out during heating; production must let it move or the RMSD is meaningless.
#
# Usage:  source /opt/md/activate.sh && bash 24_run_production.sh [ns] [concurrency]
set -uo pipefail
cd "$(dirname "$0")/.."
MD="$(pwd)"

NS="${1:-100}"
CONC="${2:-4}"
COMPLEXES="${COMPLEXES:-SUOX_Sanguinarine SUOX_Spiperone SUOX_Vorinostat}"
REPLICATES="${REPLICATES:-3}"
LOGDIR="$MD/production"
mkdir -p "$LOGDIR"

# production mdp: one per replicate, differing only in gen-seed
write_prod_mdp() {  # write_prod_mdp <file> <ns> <seed>
    local steps; steps=$(python3 -c "print(int($2*500000))")
    cat > "$1" <<EOF
integrator      = md
dt              = 0.002
nsteps          = $steps
tcoupl          = V-rescale
tc-grps         = System
tau-t           = 0.1
ref-t           = 300
pcoupl          = C-rescale
pcoupltype      = isotropic
tau-p           = 2.0
ref-p           = 1.0
compressibility = 4.5e-5
constraints     = h-bonds
constraint-algorithm = LINCS
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
nstxout-compressed      = 5000     ; every 10 ps
nstxout-compressed-type = xtc
nstenergy       = 5000
nstlog          = 5000
gen-vel         = no
; seed for the initial velocities was fixed at the NVT stage; recorded in the run manifest
EOF
}

echo "=== production: ${NS} ns x ${REPLICATES} replicates x $(echo $COMPLEXES | wc -w) complexes"
echo "    concurrency: ${CONC}"
echo "    started: $(date -Iseconds)"

MANIFEST="$LOGDIR/run_manifest.txt"
: > "$MANIFEST"
echo "# complex replicate seed ns tpr started" >> "$MANIFEST"

job=0
for cx in $COMPLEXES; do
    dir="$MD/complexes/$cx"
    [ -s "$dir/npt.gro" ] || { echo "  skip $cx (not equilibrated)"; continue; }
    for rep in $(seq 1 "$REPLICATES"); do
        seed=$((1000 + rep))
        wd="$dir/prod_rep${rep}"
        mkdir -p "$wd"
        cp "$dir/system.top" "$dir/npt.gro" "$dir/npt.cpt" "$wd/" 2>/dev/null
        write_prod_mdp "$wd/prod.mdp" "$NS" "$seed"
        # distinct initial velocities per replicate, set explicitly
        sed -i "s/gen-vel         = no/gen-vel         = yes\n gen-seed        = ${seed}/" "$wd/prod.mdp"
        (
            cd "$wd" || exit 1
            gmx grompp -f prod.mdp -c npt.gro -t npt.cpt -p system.top -o prod.tpr -maxwarn 2 \
                >grompp.log 2>&1 || { echo "  $cx rep$rep: grompp FAILED"; exit 1; }
            echo "  launched $cx rep$rep (seed $seed)"
        )
        if [ -f "$wd/prod.tpr" ]; then
            echo "$cx $rep $seed $NS $wd/prod.tpr $(date -Iseconds)" >> "$MANIFEST"
        fi
    done
done

n_jobs=$(wc -l < "$MANIFEST")
echo "  $((n_jobs - 1)) runs prepared"

# launch with a concurrency cap
running=0
while read -r cx rep seed ns tpr started; do
    case "$cx" in \#*|"") continue ;; esac
    wd=$(dirname "$tpr")
    (
        cd "$wd" || exit 1
        # mdrun with checkpointing; -maxh bounds each invocation so it resumes cleanly
        gmx mdrun -deffnm prod -nb gpu -pme gpu -bonded gpu -update gpu \
            -ntmpi 1 -ntomp "$(( $(nproc) / CONC ))" -notunepme -maxh 16.0 \
            > mdrun.log 2>&1
        echo "FINISHED $cx rep$rep rc=$?" >> "$LOGDIR/completion.txt"
    ) &
    job=$((job + 1))
    running=$((running + 1))
    if [ "$running" -ge "$CONC" ]; then
        wait -n 2>/dev/null || wait
        running=$((running - 1))
    fi
done < "$MANIFEST"

wait
echo
echo "=== all runs returned ==="
cat "$LOGDIR/completion.txt" 2>/dev/null || echo "(no completion records)"
date -Iseconds

# quick progress report: how far did each get?
echo
echo "=== trajectory progress ==="
for cx in $COMPLEXES; do
    for rep in $(seq 1 "$REPLICATES"); do
        f="$MD/complexes/$cx/prod_rep${rep}/prod.xtc"
        if [ -s "$f" ]; then
            nf=$(gmx check -f "$f" 2>&1 | grep -m1 'Step' | grep -oE 'Step[^,]*' | head -1)
            printf '  %-22s rep%s  %s\n' "$cx" "$rep" "$(gmx check -f "$f" 2>&1 | grep -m1 -E '^Step|Last frame' | sed 's/^ *//')"
        else
            printf '  %-22s rep%s  (no trajectory)\n' "$cx" "$rep"
        fi
    done
done
