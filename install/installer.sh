#!/usr/bin/env bash
# BindFlow installer.
#
# Creates two environments and links them:
#   BindFlow       workflow, FEP analysis and (optionally) MM(PB/GB)SA
#   BindFlow-TOFF  small-molecule parameterization (TOFF: OpenFF, GAFF, optionally Espaloma) and pdbfixer
#
# BindFlow runs TOFF and pdbfixer out-of-process from BindFlow-TOFF, so the heavy
# OpenMM/OpenFF/PyTorch stack never has to be solved together with the main environment.
set -euo pipefail

MAIN_ENV="BindFlow"
TOFF_ENV="BindFlow-TOFF"
MMPBSA=""
ESPALOMA=""
FORCE=""
BINDFLOW_SPEC="bindflow"
EDITABLE=""
TOFF_SPEC="toff>=0.3.0"
BINDFLOW_REF="${BINDFLOW_REF:-main}"   # git ref used to download env files if they are not next to this script

ESPALOMA_SPEC="espaloma=0.4.0"
MMPBSA_SPECS=(ambertools "mpi4py=3.1.5" tqdm)
MDRG_SPEC="git+https://github.com/ale94mleon/MDRestraintsGenerator@54ca2484801fef4f572c78639de100bf572a2581"
GMXMMPBSA_SPEC="git+https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA.git@27929e02067bc2321286809818d778a77a872010"

usage() {
    cat <<EOF
Usage: bash installer.sh [options]

  --mmpbsa           Add MM(PB/GB)SA support (AmberTools + gmx_MMPBSA) to the main environment
  --espaloma         Add Espaloma to the parameterization environment
  --name NAME        Main environment name (default: ${MAIN_ENV})
  --toff-name NAME   Parameterization environment name (default: ${TOFF_ENV})
  --bindflow SPEC    pip spec for BindFlow (default: ${BINDFLOW_SPEC}),
                     e.g. git+https://github.com/ale94mleon/BindFlow.git
  --editable DIR     Install BindFlow in editable mode from a local clone instead
  --toff SPEC        pip spec for TOFF (default: ${TOFF_SPEC})
  --force            Remove existing environments with the same names first
  -h, --help         Show this help

Requires micromamba (recommended) or mamba >= 2.
EOF
}

log() { printf '\n\033[1m==> %s\033[0m\n' "$*"; }
die() { printf '\033[31mERROR:\033[0m %s\n' "$*" >&2; exit 1; }

while [[ $# -gt 0 ]]; do
    case "$1" in
        --mmpbsa) MMPBSA=1 ;;
        --espaloma) ESPALOMA=1 ;;
        --name) MAIN_ENV="$2"; shift ;;
        --toff-name) TOFF_ENV="$2"; shift ;;
        --bindflow) BINDFLOW_SPEC="$2"; shift ;;
        --editable) EDITABLE="$(cd "$2" && pwd)"; shift ;;
        --toff) TOFF_SPEC="$2"; shift ;;
        --force) FORCE=1 ;;
        -h|--help) usage; exit 0 ;;
        *) usage; die "Unknown option: $1" ;;
    esac
    shift
done

# ---------------------------------------------------------------------------
# Conda frontend. `micromamba shell init` defines micromamba as a shell function,
# which scripts cannot see, so prefer the real binary in $MAMBA_EXE.
# ---------------------------------------------------------------------------
MM=""
for candidate in "${MAMBA_EXE:-}" "$(command -v micromamba 2>/dev/null || true)" "$(command -v mamba 2>/dev/null || true)"; do
    if [[ -n "$candidate" && -x "$candidate" ]]; then MM="$candidate"; break; fi
done
[[ -n "$MM" ]] || die "micromamba or mamba not found. Install micromamba: https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html"
echo "Using $MM ($("$MM" --version | head -n 1))"

# Use only the channels from the env files (ignores e.g. 'defaults' in ~/.condarc)
SOLVE=(--override-channels --channel-priority strict -y)

# ---------------------------------------------------------------------------
# Environment files: next to this script, or downloaded from the BindFlow repo
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
ENV_DIR="$SCRIPT_DIR/envs"
if [[ ! -f "$ENV_DIR/bindflow.yml" || ! -f "$ENV_DIR/bindflow-toff.yml" ]]; then
    ENV_DIR="$(mktemp -d)"
    for f in bindflow.yml bindflow-toff.yml; do
        curl -fsSL "https://raw.githubusercontent.com/ale94mleon/BindFlow/${BINDFLOW_REF}/install/envs/$f" -o "$ENV_DIR/$f" \
            || die "Could not download $f (ref: ${BINDFLOW_REF})"
    done
fi

env_exists() { "$MM" env list --json | grep -q "/envs/$1\""; }
prepare_env_name() {
    if env_exists "$1"; then
        [[ -n "$FORCE" ]] || die "Environment '$1' already exists. Use --force to recreate it, or choose another name."
        "$MM" env remove -n "$1" -y >/dev/null
    fi
}
prefix_of() { "$MM" run -n "$1" python -c 'import sys; print(sys.prefix)'; }

# ---------------------------------------------------------------------------
log "1/4 Creating ${TOFF_ENV} (small-molecule parameterization)"
# ---------------------------------------------------------------------------
prepare_env_name "$TOFF_ENV"
toff_extra=()
if [[ -n "$ESPALOMA" ]]; then
    toff_extra+=("$ESPALOMA_SPEC")
    # patchelf is needed below to fix Espaloma's PyTorch on Linux (see "glibc >= 2.41")
    if [[ "$(uname -s)" == "Linux" ]]; then toff_extra+=(patchelf); fi
fi
"$MM" create -n "$TOFF_ENV" -f "$ENV_DIR/bindflow-toff.yml" ${toff_extra[@]+"${toff_extra[@]}"} "${SOLVE[@]}"
# --no-deps: TOFF's own deps (rdkit, pyyaml) come from conda; pip must not replace them
"$MM" run -n "$TOFF_ENV" python -m pip install --no-deps "$TOFF_SPEC"

# glibc >= 2.41 (e.g. Ubuntu 26.04, Debian 13) refuses to load shared libraries that ask for
# an executable stack. Espaloma is limited to PyTorch 2.3 (through dgl), whose libtorch_cpu.so
# asks for one by mistake, so `import espaloma` would fail there. Clearing that flag in the
# library header fixes it; the library never runs code from the stack, so nothing else changes.
if [[ -n "$ESPALOMA" && "$(uname -s)" == "Linux" ]]; then
    TORCH_LIB="$(prefix_of "$TOFF_ENV")/lib/libtorch_cpu.so"
    if [[ -f "$TORCH_LIB" ]]; then
        "$MM" run -n "$TOFF_ENV" patchelf --clear-execstack "$TORCH_LIB"
        echo "Executable-stack flag cleared on $TORCH_LIB ($("$MM" run -n "$TOFF_ENV" patchelf --print-execstack "$TORCH_LIB"))"
    fi
fi

# ---------------------------------------------------------------------------
log "2/4 Creating ${MAIN_ENV}"
# ---------------------------------------------------------------------------
prepare_env_name "$MAIN_ENV"
main_extra=()
[[ -n "$MMPBSA" ]] && main_extra+=("${MMPBSA_SPECS[@]}")
"$MM" create -n "$MAIN_ENV" -f "$ENV_DIR/bindflow.yml" ${main_extra[@]+"${main_extra[@]}"} "${SOLVE[@]}"

pip_main() { "$MM" run -n "$MAIN_ENV" python -m pip install --no-deps "$@"; }
pip_main "$MDRG_SPEC"
[[ -n "$MMPBSA" ]] && pip_main "$GMXMMPBSA_SPEC"
if [[ -n "$EDITABLE" ]]; then pip_main -e "$EDITABLE"; else pip_main "$BINDFLOW_SPEC"; fi

# ---------------------------------------------------------------------------
log "3/4 Linking ${MAIN_ENV} -> ${TOFF_ENV}"
# ---------------------------------------------------------------------------
MAIN_PREFIX="$(prefix_of "$MAIN_ENV")"
TOFF_PREFIX="$(prefix_of "$TOFF_ENV")"
mkdir -p "$MAIN_PREFIX/etc/conda/activate.d" "$MAIN_PREFIX/etc/conda/deactivate.d"
printf 'export BINDFLOW_TOFF_PREFIX="%s"\n' "$TOFF_PREFIX" > "$MAIN_PREFIX/etc/conda/activate.d/bindflow-toff.sh"
printf 'unset BINDFLOW_TOFF_PREFIX\n' > "$MAIN_PREFIX/etc/conda/deactivate.d/bindflow-toff.sh"
echo "BINDFLOW_TOFF_PREFIX=$TOFF_PREFIX is set whenever ${MAIN_ENV} is activated"

# ---------------------------------------------------------------------------
log "4/4 Smoke test: parameterizing ethane from ${MAIN_ENV} through ${TOFF_ENV}"
# ---------------------------------------------------------------------------
FORCE_FIELDS="openff gaff${ESPALOMA:+ espaloma}"
CHECK_MMPBSA="${MMPBSA}"
export FORCE_FIELDS CHECK_MMPBSA
"$MM" run -n "$MAIN_ENV" python - <<'PYEOF'
import os
import tempfile
from pathlib import Path

import alchemlyb, MDAnalysis, MDRestraintsGenerator, bindflow  # noqa: F401
from bindflow.utils import toff_env

if os.environ.get("CHECK_MMPBSA"):
    import GMXMMPBSA.API  # noqa: F401
    print("gmx_MMPBSA: import OK")

prefix = toff_env.find_toff_prefix()
toff_env.check_toff_installation(prefix)
with tempfile.TemporaryDirectory() as tmp:
    smi = Path(tmp) / "ethane.smi"
    smi.write_text("CC\n")
    for ff in os.environ["FORCE_FIELDS"].split():
        toff_env.parameterize(prefix, smi, tmp, mol_resi_name=ff[:3].upper(), force_field_type=ff, hmr_factor=2.5)
        print(f"{ff}: OK")
print(f"BindFlow {bindflow.__version__}: installation OK")
PYEOF

log "Done"
cat <<EOF
Activate BindFlow with:

    micromamba activate ${MAIN_ENV}

You never need to activate ${TOFF_ENV}; BindFlow calls it when parameterizing small molecules.
GROMACS is not installed by this script; see the BindFlow installation docs.
EOF