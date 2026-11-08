#!/usr/bin/env bash
# 29_holo_run.sh — (re)launch the holo production grid with explicit OpenMP thread count.
#
# Why this exists
# ---------------
# The first holo launch let each mdrun pick its own thread count and GROMACS chose 16 threads.
# With nine runs at once that is 144 OpenMP threads on 32 cores: heavy oversubscription, CPU
# spinning instead of waiting on the GPU, GPU utilisation stuck around 45 % and measured throughput
# of roughly 315 ns/day. The earlier 60 ns production run used 6 threads per run and reached
# 958 ns/day on the same machine, so -ntomp is set explicitly here.
#
# All nine checkpoints exist, so every run resumes with -cpi rather than restarting.
#
# Usage: bash 29_holo_run.sh <ns> <concurrency> <ntomp>
MD="/root/olp/analysis/md"
source /opt/md/activate.sh >/dev/null 2>&1
set -uo pipefail

NS="${1:-200}"
CONC="${2:-5}"
NTOMP="${3:-6}"
COMPLEXES="${COMPLEXES:-SUOX_Sanguinarine SUOX_Spiperone SUOX_Vorinostat}"
REPLICATES="${REPLICATES:-3}"

echo "=== holo grid: ${NS} ns x ${REPLICATES} reps x 3 complexes, concurrency ${CONC}, ntomp ${NTOMP}"
echo "    started $(date -Is)"

# per-run launcher, rewritten with the chosen thread count
cat > /opt/md/holo_run_one.sh <<RUNEOF
#!/usr/bin/env bash
# one mdrun worker: resumes from prod.cpt and keeps going until the target is reached
source /opt/md/activate.sh >/dev/null 2>&1
set -uo pipefail
cd "\$1"
TARGET_PS=\$2
NTOMP=\$3

# -maxh 12 makes mdrun checkpoint and exit cleanly rather than being killed; the loop below
# re-invokes it, so a long run survives job or machine interruptions.
for attempt in \$(seq 1 60); do
    if [ -s prod.cpt ]; then
        gmx mdrun -deffnm prod -nb gpu -pme gpu -bonded gpu -update gpu \\
                  -ntomp \$NTOMP -maxh 12 -cpi prod.cpt -append >>mdrun.log 2>&1
    else
        gmx mdrun -deffnm prod -nb gpu -pme gpu -bonded gpu -update gpu \\
                  -ntomp \$NTOMP -maxh 12 >>mdrun.log 2>&1
    fi
    # last completed step in the log -> picoseconds done
    done_ps=\$(grep -E '^[[:space:]]+[0-9]+[[:space:]]+[0-9]+\\.[0-9]+' mdrun.log 2>/dev/null \\
              | tail -1 | awk '{print \$2}')
    [ -n "\$done_ps" ] || done_ps=0
    if [ "\${done_ps%.*}" -ge "\$TARGET_PS" ] 2>/dev/null; then break; fi
done
RUNEOF
chmod +x /opt/md/holo_run_one.sh

launch() {   # dir target_ps ntomp
    ( cd "$1" && setsid nohup /opt/md/holo_run_one.sh "$1" "$2" "$3" >/dev/null 2>&1 & )
}

jobs=0
for cx in $COMPLEXES; do
    for rep in 1 2 3; do
        rd="$MD/complexes/$cx/holo_rep$rep"
        [ -s "$rd/prod.tpr" ] || { echo "  skip $cx rep$rep (no tpr)"; continue; }
        # wait until the running set drops below the target concurrency
        while [ "$(pgrep -c -x gmx 2>/dev/null | head -1)" -ge "$CONC" ] 2>/dev/null; do sleep 20; done
        launch "$rd" $(( NS * 1000 )) "$NTOMP"
        echo "  launched $cx rep$rep"
        jobs=$((jobs + 1))
    done
done
echo "  ${jobs} runs launched"
