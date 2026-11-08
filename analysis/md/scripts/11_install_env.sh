#!/usr/bin/env bash
# 11_install_env.sh — install the whole MD toolchain on the rented box, then verify it.
#
# One command, idempotent, does not touch the system python.
#
#   bash 11_install_env.sh                          # default: SJTU mirror, env named 'md'
#   MIRROR=official bash 11_install_env.sh          # use the official channels
#   ENV_NAME=md2 bash 11_install_env.sh
#
# Two decisions are baked in because they were measured on the actual rented box:
#
#  1. **Mirror.** conda.anaconda.org sustained 5.4 MB/s there while the SJTU mirror sustained
#     55.9 MB/s (20 s window on the linux-64 repodata). At roughly 3 GB of packages that is the
#     difference between a ~10 minute and a ~90 minute install, which matters on a 24-hour rental.
#  2. **Use the conda that is already installed.** The provider image ships conda 24.9.2 at
#     /opt/anaconda3, and conda >= 23.10 solves with libmamba, so there is no reason to download
#     micromamba. It is only fetched as a fallback when no conda exists.
#
# The critical thing to check at the end: GROMACS must report CUDA. conda-forge's plain `gromacs`
# is a CPU build, so the script tries the CUDA variants in order and says exactly which worked.
set -uo pipefail

MD_ROOT="${MD_ROOT:-/opt/md}"
ENV_NAME="${ENV_NAME:-md}"
MIRROR="${MIRROR:-https://mirror.sjtu.edu.cn/anaconda/cloud/conda-forge}"

mkdir -p "${MD_ROOT}/bin"
CONDARC="${MD_ROOT}/condarc"
if [ "$MIRROR" = "official" ]; then
    : > "$CONDARC"
    echo "channels: official defaults"
else
    cat > "$CONDARC" <<EOF
channels:
  - $MIRROR
  - defaults
channel_priority: strict
show_channel_urls: true
EOF
    echo "channels: $MIRROR  (config: $CONDARC)"
fi
export CONDARC="$CONDARC"

# ---------------------------------------------------------------- pick a solver
if [ -x /opt/anaconda3/bin/conda ]; then
    CONDA=/opt/anaconda3/bin/conda
elif command -v conda >/dev/null 2>&1; then
    CONDA="$(command -v conda)"
else
    CONDA=""
fi
MAMBA="${MD_ROOT}/bin/micromamba"

if [ -n "$CONDA" ]; then
    echo "solver: $CONDA"
    "$CONDA" --version
    "$CONDA" config --show solver 2>/dev/null | sed 's/^/  /' || true
else
    echo "no conda found; fetching micromamba"
    arch="$(uname -m)"
    case "$arch" in
        x86_64|amd64) url="https://micro.mamba.pm/api/micromamba/linux-64/latest" ;;
        aarch64|arm64) url="https://micro.mamba.pm/api/micromamba/linux-aarch64/latest" ;;
        *) echo "unsupported arch $arch"; exit 1 ;;
    esac
    if ! curl -sSL "$url" | tar -xj -C "${MD_ROOT}" bin/micromamba; then
        echo "!! micromamba download failed - check network"; exit 1
    fi
    "$MAMBA" --version
fi

# amber tools give tleap/antechamber/sqm/parmchk2/MCPB.py (the metal site);
# gmx_mmpbsa needs a working GROMACS; mdanalysis does the trajectory analysis
COMMON_PKGS="python=3.11 ambertools gmx_mmpbsa mdanalysis parmed openbabel numpy pandas scipy matplotlib"

env_prefix() {
    if [ -n "$CONDA" ]; then
        echo "$("$CONDA" info --base 2>/dev/null)/envs/${ENV_NAME}"
    else
        echo "${MD_ROOT}/mamba/envs/${ENV_NAME}"
    fi
}

create_env() {  # create_env <gromacs spec or empty>
    local spec="$1"
    if [ -n "$CONDA" ]; then
        # shellcheck disable=SC2086
        "$CONDA" create -y -n "$ENV_NAME" $COMMON_PKGS $spec
    else
        # shellcheck disable=SC2086
        "$MAMBA" create -y -n "$ENV_NAME" $COMMON_PKGS ${spec:-gromacs}
    fi
}

remove_env() {
    if [ -n "$CONDA" ]; then "$CONDA" env remove -y -n "$ENV_NAME" >/dev/null 2>&1
    else "$MAMBA" env remove -y -n "$ENV_NAME" >/dev/null 2>&1; fi
}

# ---------------------------------------------------------------- try GROMACS variants
# conda-forge has shipped the CUDA builds under different build strings over time, so try the
# current spellings before falling back to the CPU build (kept only so the rest of the pipeline
# can still be validated).
cuda_ok=0
# NOTE the glob form: the build string is `nompi_cuda_h...`, so the pattern must be
# `*cuda*`. A bare `cuda*` matches nothing and was a real bug caught by searching the mirror
# before installing.
for entry in "nompi_cuda:gromacs=*=nompi_cuda*" "anycuda:gromacs=*=*cuda*" "cpu:gromacs"; do
    label="${entry%%:*}"; spec="${entry#*:}"
    log="${MD_ROOT}/install_${label}.log"
    echo
    echo "--- variant '$label'  ($spec)"
    if create_env "$spec" >"$log" 2>&1; then
        PREFIX="$(env_prefix)"
        gmx="${PREFIX}/bin/gmx"
        if [ -x "$gmx" ] && "$gmx" --version 2>&1 | grep -qi cuda; then
            echo "    CUDA GROMACS installed at $PREFIX"
            cuda_ok=1
            break
        fi
        if [ -x "$gmx" ]; then
            echo "    installed but CPU-only"
            [ "$label" = "cpu" ] && break
        else
            echo "    solved but no gmx binary"
        fi
        remove_env
    else
        echo "    solve/install failed (log: $log)"
        tail -5 "$log" 2>/dev/null | sed 's/^/      /'
    fi
done

PREFIX="$(env_prefix)"

# ---------------------------------------------------------------- activation helper
cat > "${MD_ROOT}/activate.sh" <<EOF
# source this to enter the MD environment
if [ -n "$CONDA" ]; then
    eval "\$("$CONDA" shell.bash hook)"
    conda activate "$PREFIX"
else
    export MAMBA_ROOT_PREFIX="${MD_ROOT}/mamba"
    eval "\$("$MAMBA" shell hook -s bash)"
    micromamba activate ${ENV_NAME}
fi
export AMBERHOME="$PREFIX"
EOF

echo
echo "=== verification (prefix $PREFIX) ==="
export PATH="${PREFIX}/bin:$PATH"

fail=0
for t in tleap antechamber sqm parmchk2 cpptraj MCPB.py obabel gmx gmx_MMPBSA; do
    if command -v "$t" >/dev/null 2>&1; then
        printf '  OK   %-12s\n' "$t"
    else
        printf '  MISS %-12s\n' "$t"; fail=$((fail+1))
    fi
done

echo
"${PREFIX}/bin/gmx" --version 2>&1 | grep -E 'GROMACS version|CUDA|GPU|SIMD' | sed 's/^/  /'
if "${PREFIX}/bin/gmx" --version 2>&1 | grep -qi cuda; then
    echo "  => GPU MD available"
else
    echo "  => !! CPU-only GROMACS; the production grid needs a CUDA build"
    echo "     try: conda install -n $ENV_NAME 'gromacs=*=cuda*'"
    echo "     or build GROMACS from source against the driver's CUDA (driver 560.35.03 -> CUDA 12.6)"
    cuda_ok=0
fi

echo
echo "AMBERHOME=${AMBERHOME:-<unset>}"
python - <<'PY'
import glob, os, sys
for r in [os.environ.get("AMBERHOME", ""), sys.prefix]:
    if r and glob.glob(os.path.join(r, "**", "leap", "cmd", "leaprc.protein.ff19SB"), recursive=True):
        print("  ff19SB data found under", r); break
else:
    print("  !! ff19SB leaprc not found - tleap will not know ff19SB")
for m in ("parmed", "MDAnalysis", "numpy"):
    try:
        mod = __import__(m); print("  ", m, getattr(mod, "__version__", "?"))
    except Exception as e:
        print("  ", m, "IMPORT FAILED:", e)
PY

echo
if [ "$cuda_ok" = 1 ] && [ "$fail" -eq 0 ]; then
    echo "ENVIRONMENT READY: CUDA GROMACS + full AmberTools."
    echo "next: source ${MD_ROOT}/activate.sh && bash 12_benchmark.sh"
else
    echo "ENVIRONMENT INCOMPLETE: missing=$fail cuda=$cuda_ok"
fi
exit 0
