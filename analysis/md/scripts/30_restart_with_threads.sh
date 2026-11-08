#!/usr/bin/env bash
# 30_restart_with_threads.sh — stop the oversubscribed grid and relaunch with a fixed thread count.
#
# The kill is done here, inside a script, rather than as an inline `pkill -f ...` over ssh: a
# `pkill -f` pattern that appears in the ssh command line also matches the remote shell running that
# very command, which terminates the session. A script's own command line is just
# "bash .../30_restart_with_threads.sh", so it can never match the pattern it is killing.
#
# The wrapper loops (./run.sh, holo_run_one.sh) must die first. They exist to restart mdrun after
# -maxh checkpoints, so killing mdrun alone just makes the wrapper start it again immediately.
MD="/root/olp/analysis/md"
source /opt/md/activate.sh >/dev/null 2>&1
set -uo pipefail

NTOMP="${1:-6}"
CONC="${2:-9}"

echo "=== stopping the old grid ==="
for pat in 'holo_run_one.sh' 'bash ./run.sh' 'run.sh 200000'; do
    pkill -f "$pat" 2>/dev/null && echo "  killed: $pat"
done
sleep 3
pkill -x gmx 2>/dev/null && echo "  killed: gmx"
sleep 8
echo "  remaining gmx: $(pgrep -c -x gmx 2>/dev/null | head -1)"
echo "  remaining wrappers: $(pgrep -fc 'holo_run_one.sh|run.sh 200000' 2>/dev/null | head -1)"

# confirm the checkpoint step so far, so the resume point is recorded
echo
echo "=== resume points (ps completed) ==="
: > /opt/md/resume_points.tsv
for d in $MD/complexes/*/holo_rep*; do
    ps=$(grep -oE 'continuing from step [0-9]+, +[0-9.]+ ps' "$d/mdrun.log" 2>/dev/null | tail -1 \
         | grep -oE '[0-9.]+ ps' | awk '{print $1}')
    [ -n "$ps" ] || ps=0
    # also take the highest step line seen, if the log has moved past the restart banner
    hi=$(grep -oE '^ *[0-9]+ +[0-9]+\.[0-9]+' "$d/mdrun.log" 2>/dev/null | tail -1 | awk '{print $2}')
    best=$(printf '%s\n%s\n' "$ps" "${hi:-0}" | sort -g | tail -1)
    echo -e "$d\t$best" >> /opt/md/resume_points.tsv
done
column -t /opt/md/resume_points.tsv | sed 's/^/  /'

echo
echo "=== relaunching with ntomp=${NTOMP}, concurrency=${CONC} ==="
cd "$MD"
COMPLEXES="SUOX_Sanguinarine SUOX_Spiperone SUOX_Vorinostat" REPLICATES=3 \
bash scripts/29_holo_run.sh 200 "$CONC" "$NTOMP"
echo
date +%s > /opt/md/base_time.txt
echo "wall-clock baseline: $(cat /opt/md/base_time.txt)"
sleep 45
echo
echo "=== after 45 s ==="
echo "  gmx processes: $(pgrep -c -x gmx 2>/dev/null | head -1)"
grep -m1 'OpenMP threads' $MD/complexes/SUOX_Sanguinarine/holo_rep1/mdrun.log | tail -1
nvidia-smi --query-gpu=utilization.gpu,power.draw --format=csv,noheader
