#!/usr/bin/env bash
# 12_benchmark.sh — measure the real ns/day on this box, at both our system sizes.
#
# Why this is the most valuable thing to do in a 24-hour rental
# ------------------------------------------------------------
# Every planning number (how long the production grid takes, how many days to rent) rests on one
# figure: ns/day on this hardware for a system of our size. Published benchmarks are only a proxy,
# so measure it. This deliberately uses pdb2gmx + a stock force field rather than the AmberTools
# route: the purpose here is to time the *hardware*, and a GROMACS-internal build cannot fail on
# ligand or metal parameterisation, so the number arrives in ~1 hour instead of a day.
#
# Force field choice barely affects throughput (PME dominates), so amber99sb-ildn is used as a
# stand-in. It is NOT the production force field (that is ff19SB via tleap) - only the timings
# carry over.
#
# Usage:
#   bash 12_benchmark.sh                 # both sizes, 2 ns each
#   NS=5 bash 12_benchmark.sh            # longer, more reliable timing
#   CONCURRENCY=4 bash 12_benchmark.sh   # also test 4 concurrent runs (the deployment plan)
set -uo pipefail

cd "$(dirname "$0")/.."
NS="${NS:-2}"
CONCURRENCY="${CONCURRENCY:-1}"
BENCH="benchmark"
mkdir -p "$BENCH"

# receptor PDBs uploaded with the project: GDF11 (901 heavy atoms) and SUOX (2895), which bracket
# our real system sizes
GDF11_PDB="${GDF11_PDB:-complexes/GDF11_GW8510/receptor.pdb}"
SUOX_PDB="${SUOX_PDB:-complexes/SUOX_Sanguinarine/receptor.pdb}"

command -v gmx >/dev/null 2>&1 || { echo "gmx not on PATH - source /opt/md/activate.sh first"; exit 1; }

echo "=== hardware ==="
nvidia-smi --query-gpu=index,name,memory.total --format=csv,noheader | sed 's/^/  /'
echo "  cpu cores: $(nproc)"

# ---------------------------------------------------------------------------- shared mdp files
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
# equilibration kept short: this is a timing exercise, not a scientific run
cat > "$BENCH/eq.mdp" <<'EOF'
integrator      = md
dt              = 0.002
nsteps          = 25000            ; 50 ps
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

build_system() {  # build_system <outdir> <pdb> <box-nm>
    local out="$1" pdb="$2" boxnm="$3"
    rm -rf "$out"; mkdir -p "$out"
    (
        cd "$out" || exit 1
        # -ignh: the uploaded receptors carry docking hydrogens; let GROMACS place its own
        gmx pdb2gmx -f "../../$pdb" -o prot.gro -p topol.top -ff amber99sb-ildn -water tip3p -ignh \
            >pdb2gmx.log 2>&1 || { echo "    pdb2gmx FAILED"; exit 1; }
        gmx editconf -f prot.gro -o box.gro -c -d "$boxnm" -bt octahedron >>build.log 2>&1
        gmx solvate -cp box.gro -cs spc216.gro -p topol.top -o solv.gro >>build.log 2>&1
        gmx grompp -f ../em.mdp -c solv.gro -p topol.top -o em.tpr -maxwarn 2 >grompp_em.log 2>&1 \
            || { echo "    grompp(em) FAILED"; exit 1; }
        gmx mdrun -deffnm em -nb cpu >mdrun_em.log 2>&1 || { echo "    EM FAILED"; exit 1; }
        gmx grompp -f ../eq.mdp -c em.gro -p topol.top -o eq.tpr -maxwarn 2 >grompp_eq.log 2>&1
        gmx mdrun -deffnm eq -nb gpu >mdrun_eq.log 2>&1 || gmx mdrun -deffnm eq >mdrun_eq.log 2>&1
        if [ -f eq.gro ]; then echo "    system built: $(grep -c . <<< "$(gmx check -f eq.gro 2>&1)") ; see atoms below"; fi
    )
    # atom count straight from the topology
    local n
    n=$(awk '/^\[ atoms \]/{f=1;next} /^\[/{f=0} f&&NF>=3{n++} END{print n+0}' "$out/topol.top" 2>/dev/null)
    echo "$n"
}

time_run() {  # time_run <dir> <ns> <label> [cpu|gpu]
    local dir="$1" ns="$2" label="$3" nb="${4:-gpu}"
    local steps; steps=$(python3 -c "print(int($ns*500000))")
    cat > "$dir/prod.mdp" <<EOF
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
    (
        cd "$dir" || exit 1
        gmx grompp -f prod.mdp -c eq.gro -p topol.top -o prod.tpr -maxwarn 2 >grompp_prod.log 2>&1 || exit 1
        local extra="-nb $nb"
        [ "$nb" = gpu ] && extra="-nb gpu -pme gpu -bonded gpu -update gpu"
        /usr/bin/time -f "%e" gmx mdrun -deffnm prod $extra -ntmpi 1 -ntomp "$(nproc)" -notunepme \
            >mdrun_prod.log 2>time.txt
        local secs; secs=$(tail -1 time.txt 2>/dev/null)
        local perf; perf=$(grep -oE '[0-9.]+ ns/day' mdrun_prod.log | tail -1)
        echo "  ${label}: ${ns} ns in ${secs}s -> ${perf:-unknown}"
    )
}

# ---------------------------------------------------------------------------- run
for spec in "gdf11|$GDF11_PDB|1.0|GDF11 (~30k atoms, no metal)" \
            "suox|$SUOX_PDB|1.0|SUOX (~55k atoms, metal site stripped for timing only)"; do
    IFS='|' read -r key pdb margin label <<< "$spec"
    echo
    echo "=== $label"
    if [ ! -f "$pdb" ]; then echo "  missing $pdb - upload the complexes/ directory first"; continue; fi
    n=$(build_system "$BENCH/$key" "$pdb" "$margin")
    echo "  atoms (solute+solvent): ${n:-?}"
    time_run "$BENCH/$key" "$NS" "$key single run" gpu
done

# ---------------------------------------------------------------------------- concurrency
if [ "$CONCURRENCY" -gt 1 ]; then
    echo
    echo "=== concurrency test: $CONCURRENCY simultaneous runs at GDF11 size"
    echo "    (the production plan runs many independent replicas; one GPU is often not saturated"
    echo "     by a 30-55k atom system, so this shows the real aggregate throughput)"
    for i in $(seq 1 "$CONCURRENCY"); do
        d="$BENCH/conc_$i"; mkdir -p "$d"
        cp "$BENCH/gdf11"/{topol.top,prod.mdp,eq.gro} "$d/" 2>/dev/null
        ( cd "$d" && gmx grompp -f prod.mdp -c eq.gro -p topol.top -o prod.tpr -maxwarn 2 \
            >grompp.log 2>&1 ) &
    done
    wait
    start=$(date +%s)
    for i in $(seq 1 "$CONCURRENCY"); do
        ( cd "$BENCH/conc_$i" && gmx mdrun -deffnm prod -nb gpu -pme gpu -bonded gpu -update gpu \
            -ntmpi 1 -ntomp "$(( $(nproc) / CONCURRENCY ))" -notunepme >mdrun.log 2>&1 ) &
    done
    wait
    end=$(date +%s)
    el=$(( end - start ))
    tot=$(python3 -c "print(round($CONCURRENCY*$NS/($el/86400.0),1))")
    single=$(grep -oE '[0-9.]+ ns/day' "$BENCH/gdf11/mdrun_prod.log" 2>/dev/null | tail -1)
    echo "  $CONCURRENCY concurrent x ${NS} ns in ${el}s -> aggregate ${tot} ns/day"
    echo "  single-run throughput was ${single:-unknown}; speed-up factor $(python3 -c "
s='${single:-0}'.split()[0]
try: print(round(($tot)/float(s),2))
except Exception: print('n/a')")"
fi

cat <<'REPORT'

==============================================================================
How to read this
==============================================================================
* The ns/day printed for each size is the number to plan with. Production needs
  4 complexes x 3 replicates x 100 ns = 1200 ns (or 2400 ns at 200 ns/replicate).
      days = 1200 / (aggregate ns/day) / 0.85        (0.85 = restarts/failures)
* SUOX is the slower size; plan on it, not on GDF11.
* If the concurrency test shows a speed-up well above 1, run one replica per GPU slot and
  raise the number of simultaneous runs rather than buying a bigger single GPU.
* Trajectory output is switched off here. Production must enable nstxout-compressed and use
  .xtc, never .trr.
REPORT
