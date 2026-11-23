#!/usr/bin/env bash
# BindFlow installer.
#
# Run it straight from GitHub (it asks about the optional parts; -y skips the questions):
#   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/ale94mleon/BindFlow/main/install/installer.sh)"
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
# The whole script is wrapped in main(), which is only called on the last line. When the script is
# downloaded (bash -c "$(curl ...)" or curl ... | bash) and the transfer is cut off, bash then never
# runs a partial script: without the final "main" line nothing happens at all.
main() {
set -euo pipefail

MAIN_ENV="BindFlow"
TOFF_ENV="BindFlow-TOFF"
MMPBSA=""
ESPALOMA=""
FORCE=""
ASSUME_YES=""
VERBOSE=""
BINDFLOW_SPEC="bindflow"
EDITABLE=""
TOFF_SPEC="toff>=0.3.0"
BINDFLOW_REF="${BINDFLOW_REF:-main}"   # git ref to download the env files from when there is no local clone (see "Environment files")

ESPALOMA_SPEC="espaloma=0.4.0"
MMPBSA_SPECS=("ambertools>=22.5" "mpi4py=3.1.5" tqdm)
MDRG_SPEC="git+https://github.com/ale94mleon/MDRestraintsGenerator@54ca2484801fef4f572c78639de100bf572a2581"
GMXMMPBSA_SPEC="git+https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA.git@27929e02067bc2321286809818d778a77a872010"

usage() {
    cat <<EOF
Usage: bash installer.sh [options]
       /bin/bash -c "\$(curl -fsSL https://raw.githubusercontent.com/ale94mleon/BindFlow/main/install/installer.sh)" -- [options]

In a terminal, the installer asks whether to add MM(PB/GB)SA and Espaloma (unless given as
options) and confirms before installing. With -y, or without a terminal (e.g. CI), it asks
nothing: only what the options request is installed.

  -y, --yes          Do not ask anything; install what the options request
  --mmpbsa           Add MM(PB/GB)SA support (AmberTools + gmx_MMPBSA) to the main environment
  --espaloma         Add Espaloma to the parameterization environment
  --name NAME        Main environment name (default: ${MAIN_ENV})
  --toff-name NAME   Parameterization environment name (default: ${TOFF_ENV})
  --bindflow SPEC    pip spec for BindFlow (default: ${BINDFLOW_SPEC}),
                     e.g. git+https://github.com/ale94mleon/BindFlow.git
  --editable DIR     Install BindFlow in editable mode from a local clone instead
  --toff SPEC        pip spec for TOFF (default: ${TOFF_SPEC})
  --force            Remove existing environments with the same names without asking
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
        -y|--yes) ASSUME_YES=1 ;;
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

# Only the path is chosen here; the file is created after "Proceed?", so declining leaves nothing behind
LOG=""
if [[ -n "$FANCY" ]]; then
    LOG="${TMPDIR:-/tmp}"
    LOG="${LOG%/}/bindflow-install-$(date +%Y%m%d-%H%M%S).log"
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
    if [[ -f "$LOG" ]]; then printf ' Log: %s' "$LOG"; fi
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
# Supported platforms: Linux (x86_64, ARM) and macOS Apple Silicon. Checked first, so an unsupported
# platform gets this explanation and not, e.g., "micromamba not found".
# ---------------------------------------------------------------------------
# Windows is not supported natively: AmberTools (antechamber for OpenFF/GAFF charges, and
# gmx_MMPBSA) and GROMACS have no Windows builds on conda-forge. WSL2 runs a real Linux, where
# everything works. Git Bash and MSYS2 report MINGW*/MSYS*, Cygwin reports CYGWIN*.
case "$(uname -s)" in
    MINGW*|MSYS*|CYGWIN*)
        die "Windows is not supported natively (AmberTools and GROMACS have no Windows builds on conda-forge). Use WSL2: run 'wsl --install' in PowerShell as administrator, then run this installer inside the Ubuntu terminal." ;;
esac

# macOS Intel (osx-64) is not supported: conda-forge's antechamber there randomly reads zero atoms from
# a valid SDF file (about 20-25% of calls, AmberTools 24.8 to 26.0). OpenFF (AM1-BCC charges) and GAFF
# both run antechamber, so parameterization fails intermittently. Upstream bug, reproduced in TOFF:
# https://github.com/ale94mleon/TOFF/blob/main/.github/scripts/antechamber-diagnostics.sh
# A shell running under Rosetta on Apple Silicon also reports x86_64 and would install osx-64 packages.
if [[ "$(uname -s)" == "Darwin" && "$(uname -m)" == "x86_64" ]]; then
    die "macOS Intel (x86_64) is not supported: AmberTools' antechamber is unreliable on this platform, so small-molecule parameterization (OpenFF and GAFF) fails intermittently. Use Linux or macOS Apple Silicon. On Apple Silicon, run the installer from a native arm64 shell, not under Rosetta."
fi

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

# Espaloma needs dgl, which conda-forge does not build for Linux ARM
ESPALOMA_AVAILABLE=1
if [[ "$(uname -s)" == "Linux" && "$(uname -m)" =~ ^(aarch64|arm64)$ ]]; then ESPALOMA_AVAILABLE=""; fi
if [[ -n "$ESPALOMA" && -z "$ESPALOMA_AVAILABLE" ]]; then
    die "--espaloma is not available on Linux ARM ($(uname -m)): its dependency dgl is not built for this platform on conda-forge. Run the installer without --espaloma."
fi

# Options for every "create": use only the channels listed in the env files (ignores e.g. 'defaults'
# in ~/.condarc), take each package from the first channel that has it (strict), and do not ask (-y).
SOLVE=(--override-channels --channel-priority strict -y)

# Helpers:
#   env_exists NAME   true if an environment called NAME already exists
#   prefix_of NAME    print the full path of environment NAME (e.g. ~/micromamba/envs/BindFlow)
env_exists() { "$MM" env list --json | grep -q "/envs/$1\""; }
prefix_of() { "$MM" run -n "$1" python -c 'import sys; print(sys.prefix)'; }

# ---------------------------------------------------------------------------
# Questions (interactive mode)
# ---------------------------------------------------------------------------
# The installer asks only when a person can answer: no -y, and the terminal (/dev/tty) can be opened.
# Questions and answers go through /dev/tty, not stdin/stdout, so they also work when the script
# itself comes through a pipe (curl ... | bash). In CI there is no terminal, so nothing is asked and
# only what the options request is installed.
INTERACTIVE=""
if [[ -z "$ASSUME_YES" ]] && { : < /dev/tty; } 2>/dev/null; then INTERACTIVE=1; fi

# ask_yes_no QUESTION DEFAULT(y|n) [IF_YES IF_NO]: true if the answer is yes; Enter takes the default
# (the capital letter in the hint). IF_YES/IF_NO, if given, are printed (dimmed) under the question
# with the outcome, so pressing Enter always shows what it chose.
ask_yes_no() {
    local hint="[y/N]" answer result note=""
    if [[ "$2" == "y" ]]; then hint="[Y/n]"; fi
    while true; do
        printf '%s?%s %s %s%s%s ' "$CYAN" "$RESET" "$1" "$DIM" "$hint" "$RESET" > /dev/tty
        read -r answer < /dev/tty || answer=""
        if [[ -z "$answer" ]]; then note=" (default)"; fi
        case "${answer:-$2}" in
            [Yy]|[Yy][Ee][Ss]) result=0; break ;;
            [Nn]|[Nn][Oo]) result=1; break ;;
        esac
        printf '  Please answer y or n.\n' > /dev/tty
    done
    if [[ $# -ge 4 ]]; then
        if [[ $result -eq 0 ]]; then answer="$3"; else answer="$4"; fi
        printf '  %s-> %s%s%s\n' "$DIM" "$answer" "$note" "$RESET" > /dev/tty
    fi
    return $result
}

printf '\n%sBindFlow installer%s\n\n' "$BOLD" "$RESET"

# Optional parts: asked only if no option decided them already
if [[ -n "$INTERACTIVE" ]]; then
    if [[ -z "$MMPBSA" ]] && ask_yes_no "Add MM(PB/GB)SA support (AmberTools + gmx_MMPBSA)?" n \
        "MM(PB/GB)SA will be installed" "No MM(PB/GB)SA"; then MMPBSA=1; fi
    if [[ -z "$ESPALOMA" ]]; then
        if [[ -n "$ESPALOMA_AVAILABLE" ]]; then
            if ask_yes_no "Add the Espaloma force field (PyTorch; a larger download)?" n \
                "Espaloma will be installed" "No Espaloma"; then ESPALOMA=1; fi
        else
            printf '%sEspaloma is not available on Linux ARM (its dependency dgl is not built for it); skipped.%s\n' "$DIM" "$RESET"
        fi
    fi
fi

# Existing environments with the same names: check both before the first (long) step, so a name
# clash never stops the installer halfway. They are removed only with --force or a "yes" here.
existing=()
for env in "$TOFF_ENV" "$MAIN_ENV"; do
    if env_exists "$env"; then existing+=("$env"); fi
done
if [[ ${#existing[@]} -gt 0 && -z "$FORCE" ]]; then
    names="$(printf '%s, ' "${existing[@]}")"; names="${names%, }"
    if [[ -n "$INTERACTIVE" ]] && ask_yes_no "These environments already exist: ${names}. Remove and recreate them?" n \
        "They will be removed and recreated" "Kept"; then
        FORCE=1
    else
        die "Environment(s) already exist: ${names}. Use --force to recreate them, or choose other names with --name/--toff-name."
    fi
fi

# ---------------------------------------------------------------------------
# Summary of what will be installed
# ---------------------------------------------------------------------------
extras=""
if [[ -n "$ESPALOMA" ]]; then extras="${extras}, Espaloma"; fi
if [[ -n "$MMPBSA" ]]; then extras="${extras}, MM(PB/GB)SA"; fi
extras="${extras#, }"
if [[ -n "$EDITABLE" ]]; then bindflow_from="$EDITABLE (editable)"; else bindflow_from="$BINDFLOW_SPEC"; fi
envs_line="${MAIN_ENV}, ${TOFF_ENV}"
if [[ ${#existing[@]} -gt 0 ]]; then envs_line="${envs_line} (existing ones are recreated)"; fi
printf '\n'
printf '  %-14s %s\n' "Frontend" "$(basename "$MM") $MM_VERSION" \
                      "Environments" "$envs_line" \
                      "Extras" "${extras:-none}" \
                      "BindFlow" "$bindflow_from"
if [[ -n "$LOG" ]]; then printf '  %-14s %s\n' "Log" "$LOG"; fi
printf '\n'
if [[ -n "$INTERACTIVE" ]]; then
    ask_yes_no "Proceed?" y || { printf 'Nothing was installed.\n'; exit 0; }
    printf '\n'
fi
if [[ -n "$LOG" ]]; then : > "$LOG"; fi
if [[ -n "$FANCY" ]]; then
    printf '%sCreating the environments downloads several GB and can take 10-20 minutes.%s\n\n' "$DIM" "$RESET"
fi

# ---------------------------------------------------------------------------
# Environment files
# ---------------------------------------------------------------------------
# Taken from the local clone when there is one: next to this script when it is run as a file
# (bash install/installer.sh), or from the --editable clone. Run from a URL (bash -c "$(curl ...)"
# or curl ... | bash) the script has no file, so they are downloaded from the BindFlow repository.
ENV_DIR=""
script_file="${BASH_SOURCE[0]:-}"
if [[ ! -f "$script_file" ]]; then script_file=""; fi
for dir in ${script_file:+"$(cd "$(dirname "$script_file")" && pwd)/envs"} ${EDITABLE:+"$EDITABLE/install/envs"}; do
    if [[ -f "$dir/bindflow.yml" && -f "$dir/bindflow-toff.yml" ]]; then ENV_DIR="$dir"; break; fi
done
if [[ -z "$ENV_DIR" ]]; then
    ENV_DIR="$(mktemp -d)"
    for f in bindflow.yml bindflow-toff.yml; do
        curl -fsSL "https://raw.githubusercontent.com/ale94mleon/BindFlow/${BINDFLOW_REF}/install/envs/$f" -o "$ENV_DIR/$f" \
            || die "Could not download $f (ref: ${BINDFLOW_REF})"
    done
fi

for env in ${existing[@]+"${existing[@]}"}; do
    run_step "Removing the existing ${env} environment" "$MM" env remove -n "$env" -y
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
MM_NAME="$(basename "$MM")"   # micromamba or mamba, whichever was used
cat <<EOF
Activate it with:

    ${BOLD}${MM_NAME} activate ${MAIN_ENV}${RESET}

You never need to activate ${TOFF_ENV}; BindFlow calls it when parameterizing small molecules.

${BOLD}Next steps${RESET}
  1. Install GROMACS. It is not part of this installer:
     https://bindflow.readthedocs.io/en/latest/source/installation.html#gromacs
  2. Check that BindFlow runs simulations with it (a few minutes):
     ${MM_NAME} install -n ${MAIN_ENV} -c conda-forge pytest
     ${MM_NAME} run -n ${MAIN_ENV} pytest --pyargs bindflow.tests.test_fep${MMPBSA:+ bindflow.tests.test_mmpbsa}
EOF
if [[ -n "$LOG" ]]; then printf '%sFull log: %s%s\n' "$DIM" "$LOG" "$RESET"; fi

}

main "$@"
