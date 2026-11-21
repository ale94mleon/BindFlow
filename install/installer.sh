#!/usr/bin/env bash
# BindFlow installer.
#
# Creates two environments and links them:
#   BindFlow       workflow, FEP analysis and (optionally) MM(PB/GB)SA
#   BindFlow-TOFF  small-molecule parameterization (TOFF: OpenFF, GAFF, optionally Espaloma) and pdbfixer
#
# BindFlow runs TOFF and pdbfixer out-of-process from BindFlow-TOFF, so the heavy
# OpenMM/OpenFF/PyTorch stack never has to be solved together with the main environment.
#
# Keep this script compatible with bash 3.2, the version macOS ships.
# Stop at the first error (-e), on undefined variables (-u) and on failures inside pipes (-o pipefail),
# so a broken step never leaves a half-installed environment that looks fine.
set -euo pipefail

MAIN_ENV="BindFlow"
TOFF_ENV="BindFlow-TOFF"
MMPBSA=""
ESPALOMA=""
FORCE=""
VERBOSE=""
BINDFLOW_SPEC="bindflow"
EDITABLE=""
TOFF_SPEC="toff>=0.3.0"
BINDFLOW_REF="${BINDFLOW_REF:-main}"   # git ref used to download env files if they are not next to this script

ESPALOMA_SPEC="espaloma=0.4.0"
MMPBSA_SPECS=("ambertools>=22.5" "mpi4py=3.1.5" tqdm)
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
  --verbose          Print the full output of every step instead of the compact progress view
  -h, --help         Show this help

Requires micromamba (recommended) or mamba >= 2.
EOF
}

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
        --verbose) VERBOSE=1 ;;
        -h|--help) usage; exit 0 ;;
        *) usage; printf 'ERROR: Unknown option: %s\n' "$1" >&2; exit 1 ;;
    esac
    shift
done

# ---------------------------------------------------------------------------
# Terminal output
# ---------------------------------------------------------------------------
# Compact view (FANCY): in an interactive terminal, every step is one line with a spinner and its
# elapsed time, plus a dimmed line below with the last line of its output. When the step ends,
# both lines are replaced by a single ✔ line. The full output of every step goes to $LOG and is
# printed only when a step fails.
# Plain view: outside a terminal (CI, pipes, files) or with --verbose, the output of every step is
# printed as it comes, under a "==> step" heading.
FANCY=""
if [[ -z "$VERBOSE" && -t 1 && "${TERM:-dumb}" != "dumb" ]]; then FANCY=1; fi

if [[ -t 1 && -z "${NO_COLOR:-}" ]]; then
    BOLD=$'\033[1m' DIM=$'\033[2m' RED=$'\033[31m' GREEN=$'\033[32m' YELLOW=$'\033[33m' CYAN=$'\033[36m' RESET=$'\033[0m'
else
    BOLD="" DIM="" RED="" GREEN="" YELLOW="" CYAN="" RESET=""
fi
# Unicode symbols only if the terminal uses UTF-8; plain ASCII otherwise
case "${LC_ALL:-${LC_CTYPE:-${LANG:-}}}" in
    *[Uu][Tt][Ff]-8*|*[Uu][Tt][Ff]8*)
        SPIN=("⠋" "⠙" "⠹" "⠸" "⠼" "⠴" "⠦" "⠧" "⠇" "⠏"); OK="✔"; FAIL="✖"; BRANCH="└" ;;
    *)
        SPIN=("|" "/" "-" "\\"); OK="+"; FAIL="x"; BRANCH="'" ;;
esac
ESC=$'\033'
CLEAR_STEP=$'\r\033[K\n\033[K\033[1A\r'   # erase the spinner line and the line below it

LOG=""
if [[ -n "$FANCY" ]]; then
    LOG="${TMPDIR:-/tmp}"
    LOG="${LOG%/}/bindflow-install-$(date +%Y%m%d-%H%M%S).log"
    : > "$LOG"
fi

die() {
    printf '%sERROR:%s %s\n' "$RED" "$RESET" "$*" >&2
    exit 1
}

fmt_time() {
    if (( $1 >= 60 )); then printf '%dm%02ds' $(($1 / 60)) $(($1 % 60)); else printf '%ds' "$1"; fi
}

# Last non-empty line written to the log by the current step, without colors or carriage-return progress bars.
# The output can contain bytes that are not valid text (e.g. a progress bar cut in the middle), so the
# filters work byte by byte (LC_ALL=C) and grep treats it as text (-a); otherwise grep prints
# "binary file matches" over the spinner.
last_line() {
    tail -n +"$((STEP_START_LINE + 1))" "$LOG" 2>/dev/null | tail -c 4000 | LC_ALL=C tr '\r\t' '\n ' \
        | LC_ALL=C sed -e "s/${ESC}\[[0-9;?]*[A-Za-z]//g" | LC_ALL=C grep -a -v -e '^[[:space:]]*$' -e '^### ' \
        | tail -n 1 || true
} 2>/dev/null

# kill_tree PID: stop PID and every process it started (micromamba, pip, ...)
kill_tree() {
    local child
    for child in $(ps -A -o pid= -o ppid= | awk -v p="$1" '$2 == p {print $1}'); do kill_tree "$child"; done
    kill "$1" 2>/dev/null || true
}

STEP_PID=""
STEP_START_LINE=0
on_exit() {
    if [[ -n "$FANCY" ]]; then printf '%s[?25h' "$ESC"; fi   # show the cursor again
}
on_interrupt() {
    if [[ -n "$STEP_PID" ]]; then kill_tree "$STEP_PID"; fi
    if [[ -n "$FANCY" ]]; then printf '%s' "$CLEAR_STEP"; fi
    printf '\n%sInterrupted.%s' "$YELLOW" "$RESET"
    if [[ -n "$LOG" ]]; then printf ' Log: %s' "$LOG"; fi
    printf '\n'
    exit 130
}
trap on_exit EXIT
trap on_interrupt INT TERM

# run_step LABEL COMMAND [ARGS...]: run one installation step (a command or a shell function).
# If it fails, the installer stops with its exit code after showing the end of its output.
run_step() {
    local label="$1"; shift
    if [[ -z "$FANCY" ]]; then
        printf '\n%s==> %s%s\n' "$BOLD" "$label" "$RESET"
        "$@"
        return
    fi

    printf '\n### %s\n### %s\n' "$label" "$*" >> "$LOG"
    STEP_START_LINE=$(wc -l < "$LOG")
    local start=$SECONDS frame=0 width line status=0 took
    # stdin from /dev/null: no step may wait for an answer that would never come
    ( "$@" ) >> "$LOG" 2>&1 < /dev/null &
    STEP_PID=$!
    printf '%s[?25l' "$ESC"   # hide the cursor while the spinner runs
    while kill -0 "$STEP_PID" 2>/dev/null; do
        width=$(tput cols 2>/dev/null || echo 80)
        line="$(last_line)"
        line="${line:0:$((width > 10 ? width - 6 : 4))}"
        printf '\r%s[K%s%s%s %s %s%s%s\n%s[K   %s%s %s%s%s[1A\r' \
            "$ESC" "$CYAN" "${SPIN[frame % ${#SPIN[@]}]}" "$RESET" "$label" \
            "$DIM" "$(fmt_time $((SECONDS - start)))" "$RESET" \
            "$ESC" "$DIM" "$BRANCH" "$line" "$RESET" "$ESC"
        frame=$((frame + 1))
        sleep 0.1
    done
    wait "$STEP_PID" || status=$?
    STEP_PID=""
    took="$(fmt_time $((SECONDS - start)))"
    printf '%s' "$CLEAR_STEP"
    if [[ $status -eq 0 ]]; then
        printf '%s%s%s %s %s%s%s\n' "$GREEN" "$OK" "$RESET" "$label" "$DIM" "$took" "$RESET"
        return
    fi
    printf '%s%s %s%s %s%s%s\n\n' "$RED" "$FAIL" "$label" "$RESET" "$DIM" "$took" "$RESET"
    tail -n +"$((STEP_START_LINE + 1))" "$LOG" | tail -n 40
    printf '\n%sThis step failed with exit code %d.%s Full log: %s\n' "$RED" "$status" "$RESET" "$LOG"
    exit "$status"
}

# ---------------------------------------------------------------------------
# Conda frontend. `micromamba shell init` defines micromamba as a shell function,
# which scripts cannot see, so prefer the real binary in $MAMBA_EXE.
# ---------------------------------------------------------------------------
# Take the first usable program, in this order: $MAMBA_EXE, micromamba on PATH, mamba on PATH.
# All later commands call it as "$MM" (e.g. "$MM" create, "$MM" run).
MM=""
for candidate in "${MAMBA_EXE:-}" "$(command -v micromamba 2>/dev/null || true)" "$(command -v mamba 2>/dev/null || true)"; do
    if [[ -n "$candidate" && -x "$candidate" ]]; then MM="$candidate"; break; fi
done
[[ -n "$MM" ]] || die "micromamba or mamba not found. Install micromamba: https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html"
MM_VERSION="$("$MM" --version | head -n 1)"

# macOS Intel (osx-64) is not supported: conda-forge's antechamber there randomly reads zero atoms from
# a valid SDF file (about 20-25% of calls, AmberTools 24.8 to 26.0). OpenFF (AM1-BCC charges) and GAFF
# both run antechamber, so parameterization fails intermittently. Upstream bug, reproduced in TOFF:
# https://github.com/ale94mleon/TOFF/blob/main/.github/scripts/antechamber-diagnostics.sh
# A shell running under Rosetta on Apple Silicon also reports x86_64 and would install osx-64 packages.
if [[ "$(uname -s)" == "Darwin" && "$(uname -m)" == "x86_64" ]]; then
    die "macOS Intel (x86_64) is not supported: AmberTools' antechamber is unreliable on this platform, so small-molecule parameterization (OpenFF and GAFF) fails intermittently. Use Linux or macOS Apple Silicon. On Apple Silicon, run the installer from a native arm64 shell, not under Rosetta."
fi

# Espaloma needs dgl, which conda-forge does not build for Linux ARM; stop before the first step
if [[ -n "$ESPALOMA" && "$(uname -s)" == "Linux" && "$(uname -m)" =~ ^(aarch64|arm64)$ ]]; then
    die "--espaloma is not available on Linux ARM ($(uname -m)): its dependency dgl is not built for this platform on conda-forge. Run the installer without --espaloma."
fi

# Options for every "create": use only the channels listed in the env files (ignores e.g. 'defaults'
# in ~/.condarc), take each package from the first channel that has it (strict), and do not ask (-y).
SOLVE=(--override-channels --channel-priority strict -y)

# ---------------------------------------------------------------------------
# Summary of what will be installed
# ---------------------------------------------------------------------------
extras=""
if [[ -n "$ESPALOMA" ]]; then extras="${extras}, Espaloma"; fi
if [[ -n "$MMPBSA" ]]; then extras="${extras}, MM(PB/GB)SA"; fi
extras="${extras#, }"
if [[ -n "$EDITABLE" ]]; then bindflow_from="$EDITABLE (editable)"; else bindflow_from="$BINDFLOW_SPEC"; fi
printf '\n%sBindFlow installer%s\n\n' "$BOLD" "$RESET"
printf '  %-14s %s\n' "Frontend" "$(basename "$MM") $MM_VERSION" \
                      "Environments" "${MAIN_ENV}, ${TOFF_ENV}" \
                      "Extras" "${extras:-none}" \
                      "BindFlow" "$bindflow_from"
if [[ -n "$LOG" ]]; then printf '  %-14s %s\n' "Log" "$LOG"; fi
printf '\n'
if [[ -n "$FANCY" ]]; then
    printf '%sCreating the environments downloads several GB and can take 10-20 minutes.%s\n\n' "$DIM" "$RESET"
fi

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

# Helpers:
#   env_exists NAME   true if an environment called NAME already exists
#   prefix_of NAME    print the full path of environment NAME (e.g. ~/micromamba/envs/BindFlow)
env_exists() { "$MM" env list --json | grep -q "/envs/$1\""; }
prefix_of() { "$MM" run -n "$1" python -c 'import sys; print(sys.prefix)'; }

# Check both names before the first (long) step, so a name clash never stops the installer halfway
for env in "$TOFF_ENV" "$MAIN_ENV"; do
    if env_exists "$env"; then
        [[ -n "$FORCE" ]] || die "Environment '$env' already exists. Use --force to recreate it, or choose another name."
        run_step "Removing the existing ${env} environment" "$MM" env remove -n "$env" -y
    fi
done

# ---------------------------------------------------------------------------
# 1. BindFlow-TOFF: small-molecule parameterization
# ---------------------------------------------------------------------------
toff_extra=()
if [[ -n "$ESPALOMA" ]]; then
    toff_extra+=("$ESPALOMA_SPEC")
    # patchelf is needed below to fix Espaloma's PyTorch on Linux (see "glibc >= 2.41")
    if [[ "$(uname -s)" == "Linux" ]]; then toff_extra+=(patchelf); fi
fi
run_step "Creating ${TOFF_ENV} (small-molecule parameterization)" \
    "$MM" create -n "$TOFF_ENV" -f "$ENV_DIR/bindflow-toff.yml" ${toff_extra[@]+"${toff_extra[@]}"} "${SOLVE[@]}"
# --no-deps: TOFF's own deps (rdkit, pyyaml) come from conda; pip must not replace them
run_step "Installing TOFF" "$MM" run -n "$TOFF_ENV" python -m pip install --no-deps "$TOFF_SPEC"

# glibc >= 2.41 (e.g. Ubuntu 26.04, Debian 13) refuses to load shared libraries that ask for
# an executable stack. Espaloma is limited to PyTorch 2.3 (through dgl), whose libtorch_cpu.so
# asks for one by mistake, so `import espaloma` would fail there. Clearing that flag in the
# library header fixes it; the library never runs code from the stack, so nothing else changes.
patch_torch() {
    local torch_lib
    torch_lib="$(prefix_of "$TOFF_ENV")/lib/libtorch_cpu.so"
    if [[ -f "$torch_lib" ]]; then
        "$MM" run -n "$TOFF_ENV" patchelf --clear-execstack "$torch_lib"
        echo "Executable-stack flag cleared on $torch_lib ($("$MM" run -n "$TOFF_ENV" patchelf --print-execstack "$torch_lib"))"
    fi
}
if [[ -n "$ESPALOMA" && "$(uname -s)" == "Linux" ]]; then
    run_step "Patching PyTorch for glibc >= 2.41" patch_torch
fi

# ---------------------------------------------------------------------------
# 2. BindFlow: the main environment
# ---------------------------------------------------------------------------
main_extra=()
if [[ -n "$MMPBSA" ]]; then main_extra+=("${MMPBSA_SPECS[@]}"); fi
run_step "Creating ${MAIN_ENV}" \
    "$MM" create -n "$MAIN_ENV" -f "$ENV_DIR/bindflow.yml" ${main_extra[@]+"${main_extra[@]}"} "${SOLVE[@]}"

# Packages that are not on conda-forge are installed with pip. --no-deps: their dependencies are
# already in the env file, and pip must not replace (or add different versions of) conda packages.
pip_main() { "$MM" run -n "$MAIN_ENV" python -m pip install --no-deps "$@"; }
run_step "Installing MDRestraintsGenerator" pip_main "$MDRG_SPEC"
if [[ -n "$MMPBSA" ]]; then run_step "Installing gmx_MMPBSA" pip_main "$GMXMMPBSA_SPEC"; fi
if [[ -n "$EDITABLE" ]]; then
    run_step "Installing BindFlow (editable)" pip_main -e "$EDITABLE"
else
    run_step "Installing BindFlow" pip_main "$BINDFLOW_SPEC"
fi

# ---------------------------------------------------------------------------
# 3. Link BindFlow -> BindFlow-TOFF
# ---------------------------------------------------------------------------
# BindFlow must know where BindFlow-TOFF is. micromamba/mamba/conda run every script in an
# environment's etc/conda/activate.d/ when it is activated (and deactivate.d/ when deactivated).
# So we write a one-line script there that sets BINDFLOW_TOFF_PREFIX to the full path of
# BindFlow-TOFF; bindflow.utils.toff_env.find_toff_prefix() reads it. Without activation (e.g. running
# <env>/bin/python directly), BindFlow falls back to looking for BindFlow-TOFF next to BindFlow.
link_envs() {
    local main_prefix toff_prefix
    main_prefix="$(prefix_of "$MAIN_ENV")"
    toff_prefix="$(prefix_of "$TOFF_ENV")"
    mkdir -p "$main_prefix/etc/conda/activate.d" "$main_prefix/etc/conda/deactivate.d"
    printf 'export BINDFLOW_TOFF_PREFIX="%s"\n' "$toff_prefix" > "$main_prefix/etc/conda/activate.d/bindflow-toff.sh"
    printf 'unset BINDFLOW_TOFF_PREFIX\n' > "$main_prefix/etc/conda/deactivate.d/bindflow-toff.sh"
    echo "BINDFLOW_TOFF_PREFIX=$toff_prefix is set whenever ${MAIN_ENV} is activated"
}
run_step "Linking ${MAIN_ENV} to ${TOFF_ENV}" link_envs

# ---------------------------------------------------------------------------
# 4. Smoke test
# ---------------------------------------------------------------------------
# Run a short Python script inside BindFlow: import the main packages, find BindFlow-TOFF, check
# its TOFF version and parameterize ethane with each installed force field. Each check prints
# "<name>: OK"; the compact view lists them under the step. The settings are passed to Python as
# environment variables, and PYTHONUNBUFFERED makes the lines appear as soon as they are printed.
smoke_test() {
    FORCE_FIELDS="openff gaff${ESPALOMA:+ espaloma}" CHECK_MMPBSA="${MMPBSA}" PYTHONUNBUFFERED=1 \
        "$MM" run -n "$MAIN_ENV" python - <<'PYEOF'
import os
import tempfile
from pathlib import Path

import alchemlyb, MDAnalysis, MDRestraintsGenerator, bindflow  # noqa: F401
from bindflow.utils import toff_env

print(f"BindFlow {bindflow.__version__} import: OK")
if os.environ.get("CHECK_MMPBSA"):
    import GMXMMPBSA.API  # noqa: F401
    print("gmx_MMPBSA import: OK")

prefix = toff_env.find_toff_prefix()
version = toff_env.check_toff_installation(prefix)
print(f"TOFF {version}: OK")
with tempfile.TemporaryDirectory() as tmp:
    smi = Path(tmp) / "ethane.smi"
    smi.write_text("CC\n")
    for ff in os.environ["FORCE_FIELDS"].split():
        resname = ff[:3].upper()
        try:
            toff_env.parameterize(prefix, smi, tmp, mol_resi_name=resname, force_field_type=ff, hmr_factor=2.5)
        except RuntimeError:
            # The error only shows the end of the output; print the whole log before tmp is deleted
            log_file = Path(tmp) / f"{resname}_toff.log"
            if log_file.is_file():
                print(f"--- {log_file.name} ---\n{log_file.read_text()}")
            raise
        print(f"Parameterization with {ff}: OK")
PYEOF
}
run_step "Testing the installation" smoke_test
if [[ -n "$FANCY" ]]; then
    tail -n +"$((STEP_START_LINE + 1))" "$LOG" | sed -n 's/: OK$//p' | while IFS= read -r check; do
        printf '    %s%s%s %s\n' "$GREEN" "$OK" "$RESET" "$check"
    done
fi

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------
printf '\n%s%s BindFlow is installed%s %s(%s)%s\n\n' "$BOLD$GREEN" "$OK" "$RESET" "$DIM" "$(fmt_time "$SECONDS")" "$RESET"
cat <<EOF
Activate it with:

    ${BOLD}micromamba activate ${MAIN_ENV}${RESET}

You never need to activate ${TOFF_ENV}; BindFlow calls it when parameterizing small molecules.
GROMACS is not installed by this script; see the BindFlow installation docs.
EOF
if [[ -n "$LOG" ]]; then printf '%sFull log: %s%s\n' "$DIM" "$LOG" "$RESET"; fi
