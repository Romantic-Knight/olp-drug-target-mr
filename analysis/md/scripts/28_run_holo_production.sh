#!/usr/bin/env bash
# 28_run_holo_production.sh — production MD with the cofactor present.
#
# This is the experiment the earlier run could not answer. The first pass used a protein-only
# receptor, so a ligand that left its pose might simply have fallen into the emptied cofactor
# cavity; MTE and the catalytic sulfate are now in the pocket (see 27_build_suox_holo.sh), which
# removes that explanation. Whatever these trajectories show about pose stability is therefore a
# statement about the ligand in a pocket that contains the cofactor - still without the Mo, which
# cannot be parameterised here.
#
# 200 ns per replicate rather than 60 ns: with the box already equilibrated the marginal cost is
# small, and 60 ns was short enough that several runs had not clearly plateaued.
#
# All-atom compressed trajectories are written every 10 ps and stripped to solute afterwards, so
# the download stays manageable (the same approach used for the first pass).
#
# Usage: bash 28_run_holo_production.sh <ns> <concurrency>
MD="/root/olp/analysis/md"
source /opt/md/activate.sh >/dev/null 2>&1
set -uo pipefail

NS="${1:-200}"
CONC="${2:-5}"
COMPLEXES="${COMPLEXES:-SUOX_Sanguinarine SUOX_Spiperone SUOX_Vorinostat}"
REPLICATES="${REPLICATES:-3}"
SEEDS="${SEEDS:-2001 2002 2003}"

TOTAL_NS=$(echo "$NS * $REPLICATES * $(echo $COMPLEXES | wc -w)" | bc)
echo "=== holo production: ${NS} ns x ${REPLICATES} replicates x $(echo $COMPLEXES | wc -w) complexes = ${TOTAL_NS} ns"
echo "    concurrency ${CONC}, started $(date -Is)"

jobs=0
for cx in $COMPLEXES; do
    dir="$MD/complexes/$cx"
    [ -s "$dir/hnpt.gro" ] || { echo "  skip $cx (not equilibrated)"; continue; }
    [ -s "$dir/holo.top" ] || { echo "  skip $cx (no holo.top)"; continue; }

    for rep in $(seq 1 "$REPLICATES"); do
        seed=$(echo "$SEEDS" | cut -d' ' -f"$rep")
        rd="$dir/holo_rep$rep"
        mkdir -p "$rd"
        cp -n "$dir/holo.top" "$rd/holo.top" 2>/dev/null || true

        if [ ! -s "$rd/prod.mdp" ]; then
            steps=$(( NS * 500000 ))          # 2 fs timestep
            cat > "$rd/prod.mdp" <<EOF
define              = -DPOSRES
integrator          = md
dt                  = 0.002
nsteps              = ${steps}
nstxout-compressed  = 5000
nstenergy           = 5000
nstlog              = 50000
tcoupl              = V-rescale
tc-grps             = System
tau-t               = 0.1
ref-t               = 300
pcoupl              = C-rescale
pcoupltype          = isotropic
tau-p                 = 2.0
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
            printf 'Protein\nSystem\n' | gmx grompp -f "$rd/prod.mdp" -c "$dir/hnpt.gro" -r "$dir/hnpt.gro" \
                -t "$dir/hnpt.cpt" -p "$dir/holo.top" -o "$rd/prod.tpr" -maxwarn 2 >"$rd/grompp.log" 2>&1
            if [ ! -s "$rd/prod.tpr" ]; then
                echo "  $cx rep$rep: grompp FAILED"; grep -iE 'error|fatal' "$rd/grompp.log" | head -3 | sed 's/^/      /'
                continue
            fi
        fi

        # A per-run launcher script instead of an inline `bash -c "..."`. The inline form broke
        # because the outer shell expanded $APPEND inside the double-quoted string while `set -u`
        # was active, giving "APPEND: unbound variable"; writing the logic to a file removes the
        # nested quoting entirely.
        #
        # -maxh 12 makes mdrun checkpoint and return cleanly, and the loop re-invokes it until the
        # requested length is reached, so a long run survives a job or machine stop.
        cat > "$rd/run.sh" <<'RUNEOF'
#!/usr/bin/env bash
source /opt/md/activate.sh >/dev/null 2>&1
set -uo pipefail
cd "$(dirname "$0")"
TARGET_PS=$1
for attempt in $(seq 1 40); do
    if [ -s prod.cpt ]; then
        gmx mdrun -deffnm prod -nb gpu -pme gpu -bonded gpu -update gpu \
                  -maxh 12 -cpi prod.cpt -append >>mdrun.log 2>&1
    else
        gmx mdrun -deffnm prod -nb gpu -pme gpu -bonded gpu -update gpu \
                  -maxh 12 >>mdrun.log 2>&1
    fi
    done_ps=$(gmx check -f prod.xtc 2>&1 | awk '/Last frame/{print $4}' || echo 0)
    [ -n "$done_ps" ] || done_ps=0
    if [ "$done_ps" -ge "$TARGET_PS" ] 2>/dev/null; then break; fi
done
RUNEOF
        chmod +x "$rd/run.sh"
        ( cd "$rd" && setsid nohup ./run.sh $(( NS * 1000 )) >/dev/null 2>&1 & )
        echo "  launched $cx rep$rep (seed $seed, ${NS} ns)"
        jobs=$((jobs + 1))
        while [ "$(pgrep -c -x mdrun 2>/dev/null | head -1)" -ge "$CONC" ] 2>/dev/null; do sleep 30; done
    done
done

echo "  ${jobs} runs prepared"
