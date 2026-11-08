#!/usr/bin/env bash
# 00_setup_server.sh — sanity-check the MD environment before any real work.
#
# Run after `pixi install` (or `mamba env create -f pixi.toml`). This installs nothing; it only
# verifies that every tool the pipeline needs is present and reports its version, so a missing
# piece surfaces now rather than 12 hours into a GPU run.
set -euo pipefail

echo "=== MD environment check ==="
fail=0

check() {  # check <command> <version-flag> <label>
    local cmd="$1" flag="$2" label="$3"
    if command -v "$cmd" >/dev/null 2>&1; then
        printf '  OK   %-14s %s\n' "$label" "$("$cmd" "$flag" 2>&1 | head -1)"
    else
        printf '  MISS %-14s not on PATH\n' "$label"
        fail=$((fail + 1))
    fi
}

check tleap        -h "AmberTools tleap"
check antechamber  -h "antechamber"
check parmchk2     -h "parmchk2"
check sqm          -h "sqm"
check gmx          --version "GROMACS"
check cpptraj      --version "cpptraj"
check obabel       -V "OpenBabel"

# MCPB.py ships inside AmberTools as a script rather than a binary
if command -v MCPB.py >/dev/null 2>&1; then
    echo "  OK   MCPB.py        present"
else
    echo "  MISS MCPB.py        required for the Mo-MPT bonded model"
    fail=$((fail + 1))
fi

if command -v gmx_MMPBSA >/dev/null 2>&1; then
    echo "  OK   gmx_MMPBSA    $(gmx_MMPBSA --version 2>&1 | head -1)"
else
    echo "  MISS gmx_MMPBSA    required for MM-PBSA"
    fail=$((fail + 1))
fi

# CUDA: without it the production grid is not affordable, so flag it loudly
if command -v nvidia-smi >/dev/null 2>&1; then
    echo "  GPU:"
    nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader | sed 's/^/    /'
    if gmx --version 2>&1 | grep -qi cuda; then
        echo "    GROMACS reports CUDA support: good"
    else
        echo "    WARNING: this GROMACS build does not report CUDA. Reinstall a CUDA build"
        echo "             (e.g. gromacs=*=nompi_cuda*) before the production grid."
    fi
else
    echo "  WARNING: no nvidia-smi; GPU MD will not be available"
fi

echo
if [ "$fail" -eq 0 ]; then
    echo "All required tools present."
else
    echo "$fail tool(s) missing - fix before continuing."
fi

echo
echo "=== force field data ==="
echo "AMBERHOME=${AMBERHOME:-<unset>}"
python - <<'PY' || true
import os, glob, sys
roots = [os.environ.get("AMBERHOME", ""), sys.prefix]
for r in roots:
    if not r:
        continue
    hits = glob.glob(os.path.join(r, "**", "leap", "cmd", "leaprc.protein.ff19SB"),
                     recursive=True)
    if hits:
        print("  ff19SB found under", r)
        break
else:
    print("  ff19SB leaprc NOT found - tleap will not know ff19SB")
PY
exit "$fail"
