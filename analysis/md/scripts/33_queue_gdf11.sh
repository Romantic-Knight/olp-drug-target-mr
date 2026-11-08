#!/usr/bin/env bash
# 33_queue_gdf11.sh — start the GDF11 grid once the SUOX grid has finished.
#
# Why queue instead of launching now
# ----------------------------------
# Nine SUOX runs at -ntomp 6 already use 54 OpenMP threads on 32 cores; adding three more would
# push that to 72 and reintroduce the oversubscription that cost us a factor of 2.3 earlier
# (386 -> 907 ns/day). Waiting costs nothing here: SUOX needs about two days of the seven-day
# window and GDF11 about half a day.
#
# The wait is a poll on the number of live gmx processes, not a fixed sleep, so it adapts.
#
# Usage: bash 33_queue_gdf11.sh [ns] [concurrency] [ntomp]
MD="/root/olp/analysis/md"
source /opt/md/activate.sh >/dev/null 2>&1
set -uo pipefail

NS="${1:-200}"
CONC="${2:-3}"
NTOMP="${3:-6}"

echo "=== queued: GDF11 grid after SUOX, ${NS} ns x 3 replicates, ntomp ${NTOMP}"
date -Is

# wait for the SUOX grid to drain
while :; do
    live=$(pgrep -c -x gmx 2>/dev/null | head -1)
    [ -n "$live" ] || live=0
    if [ "$live" -eq 0 ]; then
        echo "  SUOX grid finished at $(date -Is)"
        break
    fi
    sleep 120
done

# also require the SUOX trajectories to have reached the requested length, otherwise the queue
# would fire while runs are merely between mdrun invocations
cd "$MD"
for d in complexes/SUOX_*/holo_rep*; do
    if [ -s "$d/prod.xtc" ]; then
        # count frames: 10 ps per frame with nstxout-compressed 5000 at 2 fs
        frames=$(gmx check -f "$d/prod.xtc" 2>&1 | awk '/Last frame/{print $NF}' | tr -d ',')
        echo "  $(basename $(dirname $d))/$(basename $d): ${frames:-?} frames"
    fi
done

echo
echo "=== launching GDF11 production ==="
dir="$MD/complexes/GDF11_GW8510"
[ -s "$dir/gnpt.gro" ] || { echo "  GDF11 not equilibrated - aborting"; exit 1; }

# per-run launcher, same pattern as the SUOX grid: resume from checkpoint, loop until the target
cat > "$dir/gdf11_run_one.sh" <<'RUNEOF'
#!/usr/bin/env bash
source /opt/md/activate.sh >/dev/null 2>&1
set -uo pipefail
cd "$1"
TARGET_PS=$2
NTOMP=$3
for attempt in $(seq 1 60); do
    if [ -s prod.cpt ]; then
        gmx mdrun -deffnm prod -nb gpu -pme gpu -bonded gpu -update gpu \
                  -ntomp $NTOMP -maxh 12 -cpi prod.cpt -append >>mdrun.log 2>&1
    else
        gmx mdrun -deffnm prod -nb gpu -pme gpu -bonded gpu -update gpu \
                  -ntomp $NTOMP -maxh 12 >>mdrun.log 2>&1
    fi
    done_ps=$(grep -E '^[[:space:]]+[0-9]+[[:space:]]+[0-9]+\.[0-9]+' mdrun.log 2>/dev/null | tail -1 | awk '{print $2}')
    [ -n "$done_ps" ] || done_ps=0
    if [ "${done_ps%.*}" -ge "$TARGET_PS" ] 2>/dev/null; then break; fi
done
RUNEOF
chmod +x "$dir/gdf11_run_one.sh"

for rep in 1 2 3; do
    rd="$dir/gdf11_rep$rep"
    mkdir -p "$rd"
    cp -n "$dir/gdf11.top" "$rd/gdf11.top" 2>/dev/null || true
    seed=$(( 3000 + rep ))
    if [ ! -s "$rd/prod.mdp" ]; then
        cat > "$rd/prod.mdp" <<EOF
define              = -DPOSRES
integrator          = md
dt                  = 0.002
nsteps              = $(( NS * 500000 ))
nstxout-compressed  = 5000
nstenergy           = 5000
nstlog              = 50000
tcoupl              = V-rescale
tc-grps             = System
tau-t               = 0.1
ref-t               = 300
pcoupl              = C-rescale
pcoupltype          = isotropic
tau-p               = 2.0
ref-p               = 1.0
compressibility     = 4.5e-5
constraints         = h-bonds
constraint-algorithm = LINCS
cutoff-scheme       = Verlet
coulombtype         = PME
rcoulomb            = 1.0
rvdw                = 1.0
gen-vel             = yes
gen-temp            = 300
gen-seed            = ${seed}
EOF
        printf 'Protein\nSystem\n' | gmx grompp -f "$rd/prod.mdp" -c "$dir/gnpt.gro" -r "$dir/gnpt.gro" \
            -t "$dir/gnpt.cpt" -p "$dir/gdf11.top" -o "$rd/prod.tpr" -maxwarn 2 >"$rd/grompp.log" 2>&1
        if [ ! -s "$rd/prod.tpr" ]; then
            echo "  rep$rep grompp FAILED"; grep -iE 'error|fatal' "$rd/grompp.log" | head -3 | sed 's/^/      /'
            continue
        fi
    fi
    ( cd "$rd" && setsid nohup "$dir/gdf11_run_one.sh" "$rd" $(( NS * 1000 )) "$NTOMP" >/dev/null 2>&1 & )
    echo "  launched GDF11 rep$rep (seed $seed, ${NS} ns)"
done
echo "=== GDF11 queue complete at $(date -Is) ==="
