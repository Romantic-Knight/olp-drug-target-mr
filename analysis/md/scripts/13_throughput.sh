#!/usr/bin/env bash
# 13_throughput.sh — the planning numbers: system-size scaling and GPU concurrency.
#
# 12_benchmark.sh gave the first real number (a GDF11-sized system at full GPU offload). This script
# answers the two questions that decide how the GPU hours are spent:
#
#   1. **Size scaling.** SUOX cannot go through pdb2gmx (MTE/MO/SO4 have no GROMACS parameters), so
#      its timing proxy is the same protein with hetero atoms stripped. The cofactors are ~35 atoms
#      out of ~55,000 - a fair throughput proxy. Real SUOX runs use the AmberTools route anyway.
#
#   2. **Concurrency.** We have 12 independent runs (4 complexes x 3 replicates) and ONE GPU. A
#      30-55k atom system does not saturate a 4090, so what matters is how many simultaneous mdrun
#      processes maximise *aggregate* ns/day. That is worth more than any hardware upgrade.
#
# Usage:  source /opt/md/activate.sh && bash 13_throughput.sh
set -uo pipefail
cd "$(dirname "$0")/.."          # -> analysis/md
MD="$(pwd)"
BENCH="$MD/benchmark"
GPU_FLAGS="-nb gpu -pme gpu -bonded gpu -update gpu"
NS="${NS:-1}"
mkdir -p "$BENCH"

write_mdp() {
cat > "$BENCH/em.mdp" <<'EOF'
integrator      = steep
nsteps          = 20000
emtol           = 1000.0
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
pbc             = xyz
EOF
cat > "$BENCH/eq.mdp" <<'EOF'
integrator      = md
dt              = 0.002
nsteps          = 10000
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
gen-vel         = yes
gen-temp        = 300
nstxout-compressed = 0
EOF
}

write_prod() {  # write_prod <dir> <ns>
    local steps; steps=$(python3 -c "print(int($2*500000))")
    cat > "$1/prod.mdp" <<EOF
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
nstxout-compressed = 0
nstlog          = 5000
EOF
}

build() {  # build <dir> <pdb(absolute)> <strip hetero yes/no>
    local dir="$1" pdb="$2" strip="$3"
    rm -rf "$dir"; mkdir -p "$dir"
    [ -f "$pdb" ] || { echo "MISSING_INPUT"; return 1; }
    local src="$pdb"
    if [ "$strip" = "yes" ]; then
        grep -E '^(ATOM|TER)' "$src" > "$dir/_protein_only.pdb"
        printf 'END\n' >> "$dir/_protein_only.pdb"
        src="$dir/_protein_only.pdb"
    fi
    (
      cd "$dir" || exit 1
      gmx pdb2gmx -f "$src" -o prot.gro -p topol.top -ff amber99sb-ildn -water tip3p -ignh \
          >pdb2gmx.log 2>&1 || exit 2
      gmx editconf -f prot.gro -o box.gro -c -d 1.0 -bt octahedron >build.log 2>&1 || exit 3
      gmx solvate -cp box.gro -cs spc216.gro -p topol.top -o solv.gro >>build.log 2>&1 || exit 4
      gmx grompp -f "$BENCH/em.mdp" -c solv.gro -p topol.top -o em.tpr -maxwarn 2 >grompp_em.log 2>&1 || exit 5
      gmx mdrun -deffnm em -nb cpu >mdrun_em.log 2>&1 || exit 6
      gmx grompp -f "$BENCH/eq.mdp" -c em.gro -p topol.top -o eq.tpr -maxwarn 2 >grompp_eq.log 2>&1 || exit 7
      gmx mdrun -deffnm eq $GPU_FLAGS >mdrun_eq.log 2>&1 || gmx mdrun -deffnm eq >mdrun_eq.log 2>&1 || exit 8
      gmx dump -s eq.tpr 2>/dev/null | sed -n 's/.*natoms *= *\([0-9][0-9]*\).*/\1/p' | head -1 > natoms.txt
      exit 0
    )
    local rc=$?
    if [ $rc -ne 0 ]; then
        echo "BUILD_FAILED(step $rc)"
        tail -3 "$dir/pdb2gmx.log" 2>/dev/null | sed 's/^/      /'
        return 1
    fi
    cat "$dir/natoms.txt" 2>/dev/null | tr -d '\n'
}

timed_run() {  # timed_run <dir> <ns> <tprname>  -> prints wall seconds
    local dir="$1" ns="$2" tpr="${3:-run}"
    write_prod "$dir" "$ns"
    ( cd "$dir" && gmx grompp -f prod.mdp -c eq.gro -p topol.top -o "$tpr.tpr" -maxwarn 2 \
        >"grompp_$tpr.log" 2>&1 ) || { echo ERR; return 1; }
    [ -f "$dir/$tpr.tpr" ] || { echo ERR; return 1; }
    local t0 t1
    t0=$(date +%s.%N)
    ( cd "$dir" && gmx mdrun -deffnm "$tpr" $GPU_FLAGS -ntmpi 1 -ntomp "$(nproc)" -notunepme \
        >"mdrun_$tpr.log" 2>&1 )
    t1=$(date +%s.%N)
    if ! grep -q 'Performance:' "$dir/mdrun_$tpr.log" 2>/dev/null; then
        echo "mdrun failed:" >&2; tail -3 "$dir/mdrun_$tpr.log" >&2; echo ERR; return 1
    fi
    awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.2f", b-a}'
}

nsday() { awk -v s="$1" -v n="$2" 'BEGIN{ if(s>0) printf "%.0f", n/s*86400; else print "?" }'; }

echo "=== throughput on $(nvidia-smi --query-gpu=name --format=csv,noheader | head -1) ==="
echo "  driver       : $(nvidia-smi --query-gpu=driver_version --format=csv,noheader | head -1)"
echo "  GROMACS      : $(gmx -quiet --version 2>/dev/null | grep -m1 GROMACS | sed 's/.*GROMACS //')"
echo "  gpu offload  : $GPU_FLAGS"
echo "  measurement  : ${NS} ns per run"
echo

write_mdp
SUM="$BENCH/throughput_summary.txt"
: > "$SUM"

echo "=== size scaling ==="
declare -A WALL
for spec in "gdf11|complexes/GDF11_GW8510/receptor.pdb|no|GDF11" \
            "suox|complexes/SUOX_Sanguinarine/receptor.pdb|yes|SUOX_stripped"; do
    IFS='|' read -r key relpdb strip label <<< "$spec"
    abs="$MD/$relpdb"
    echo "--- building $label  ($relpdb)"
    nat=$(build "$BENCH/$key" "$abs" "$strip")
    case "$nat" in
        MISSING_INPUT) echo "    input not found: $abs"; continue ;;
        BUILD_FAILED*|"") echo "    build failed, skipping"; continue ;;
    esac
    echo "    atoms: $nat"
    secs=$(timed_run "$BENCH/$key" "$NS" run)
    [ "$secs" = "ERR" ] && { echo "    run failed, skipping"; continue; }
    WALL[$key]="$secs"
    printf '    %-16s %8s atoms  %8.2f s  ->  %6s ns/day\n' "$label" "$nat" "$secs" "$(nsday "$secs" "$NS")"
    echo "$label atoms=$nat wall_s=$secs nsday=$(nsday "$secs" "$NS")" >> "$SUM"
done

echo
echo "=== concurrency: N simultaneous GDF11-size runs, ${NS} ns each ==="
SRC="$BENCH/gdf11"
if [ -f "$SRC/eq.gro" ]; then
    for N in 1 2 4 6; do
        for i in $(seq 1 "$N"); do
            d="$BENCH/c${N}_$i"; rm -rf "$d"; mkdir -p "$d"
            cp "$SRC"/{eq.gro,topol.top} "$d/" 2>/dev/null
            write_prod "$d" "$NS"
            ( cd "$d" && gmx grompp -f prod.mdp -c eq.gro -p topol.top -o run.tpr -maxwarn 2 \
                >grompp.log 2>&1 )
        done
        t0=$(date +%s.%N)
        for i in $(seq 1 "$N"); do
            ( cd "$BENCH/c${N}_$i" && gmx mdrun -deffnm run $GPU_FLAGS -ntmpi 1 \
                -ntomp "$(( $(nproc) / N ))" -notunepme >mdrun.log 2>&1 ) &
        done
        wait
        t1=$(date +%s.%N)
        secs=$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.2f", b-a}')
        per=$(nsday "$secs" "$NS")
        agg=$(awk -v n="$N" -v s="$secs" -v d="$NS" 'BEGIN{ if(s>0) printf "%.0f", n*d/s*86400; else print "?"}')
        printf '    N=%-2s  wall %7.2f s   per-run %6s ns/day   AGGREGATE %6s ns/day\n' "$N" "$secs" "$per" "$agg"
        echo "concurrency N=$N wall_s=$secs per_run_nsday=$per aggregate_nsday=$agg" >> "$SUM"
    done
else
    echo "    GDF11 build unavailable; concurrency sweep skipped"
fi

echo
echo "=== summary written to $SUM ==="
cat "$SUM"
