#!/usr/bin/env bash
# 10_inventory.sh — first thing to run on the rented box.
#
# Prints every hardware number that decides how the project is planned: how many GPUs, whether the
# GROMACS build will actually have CUDA, how many CPU cores the MM-PBSA stage can use, how much RAM
# the concurrent runs can occupy, and whether the disk will survive. Run this before installing
# anything: if the machine is not what was ordered, that is worth knowing in the first minute, not
# after a day of setup.
set -uo pipefail

echo "=================== rented instance inventory ==================="
date
echo "host: $(hostname)   uptime:$(uptime -p 2>/dev/null | sed 's/up //')"

echo
echo "--- GPU ---"
if command -v nvidia-smi >/dev/null 2>&1; then
    nvidia-smi --query-gpu=index,name,memory.total,driver_version,compute_cap \
               --format=csv,noheader | nl -w2 -s') '
    echo "gpu count: $(nvidia-smi --query-gpu=name --format=csv,noheader | wc -l)"
    echo
    nvidia-smi | tail -n +9 | head -n 14
else
    echo "nvidia-smi NOT FOUND - no NVIDIA driver; GPU MD will be impossible"
fi

echo
echo "--- CUDA toolkit ---"
if command -v nvcc >/dev/null 2>&1; then nvcc --version | tail -2; else echo "nvcc absent (not needed if GROMACS comes from conda)"; fi

echo
echo "--- CPU ---"
lscpu | grep -E '^(Model name|CPU\(s\)|Thread|Core|Socket|MHz|NUMA node\(s\))' | sed 's/^/  /'

echo
echo "--- memory ---"
free -h | sed 's/^/  /'

echo
echo "--- disk ---"
df -hT | grep -vE 'tmpfs|udev|loop' | sed 's/^/  /'
echo "  note: a data disk that is NOT the system disk is the one likely to survive a rebuild."

echo
echo "--- OS / kernel ---"
echo "  $(. /etc/os-release 2>/dev/null && echo "$PRETTY_NAME")  kernel $(uname -r)  arch $(uname -m)"

echo
echo "--- preinstalled MD tools (unlikely, but check before installing) ---"
for t in gmx gmx_mpi tleap antechamber sqm cpptraj parmchk2 obabel MCPB.py conda mamba micromamba pixi; do
    if command -v "$t" >/dev/null 2>&1; then
        printf '  present: %-12s %s\n' "$t" "$(command -v "$t")"
    fi
done

echo
echo "--- network reachability (conda-forge, needed for the install) ---"
if curl -sS --max-time 10 -o /dev/null -w '  conda.anaconda.org HTTP %{http_code} in %{time_total}s\n' \
        https://conda.anaconda.org/conda-forge/noarch/repodata.json 2>/dev/null; then :; else
    echo "  WARNING: cannot reach conda-forge; the install step will fail"
fi

echo
echo "=================== what to check against the order ==================="
cat <<'CHECKLIST'
  [ ] GPU model and count match what was purchased
  [ ] driver_version new enough for the GROMACS build (r525+ for CUDA 12; r570+ for Blackwell/5090)
  [ ] CPU cores >= 64        (gmx_MMPBSA over ~12,000 frames is core-parallel)
  [ ] RAM >= 256 GB          (12 concurrent mdrun + PME + analysis)
  [ ] free disk >= 1 TB on a volume that will persist
  [ ] a usable CUDA GROMACS is obtainable (next script installs it)
CHECKLIST
