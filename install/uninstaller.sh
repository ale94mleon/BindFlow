#!/usr/bin/env bash
# BindFlow uninstaller: removes the environments created by installer.sh and, if asked, the files
# that Espaloma leaves in the home directory.
#
# Run it straight from GitHub (it lists what it will remove and asks first; -y skips the question):
#   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/ale94mleon/BindFlow/main/install/uninstaller.sh)"
#
# Keep this script compatible with bash 3.2, the version macOS ships.
set -euo pipefail

MAIN_ENV="BindFlow"
TOFF_ENV="BindFlow-TOFF"
TOFF_ENV_GIVEN=""
ASSUME_YES=""
DRY_RUN=""
CACHES=""

usage() {
    cat <<EOF
Usage: bash uninstaller.sh [options]
       /bin/bash -c "\$(curl -fsSL https://raw.githubusercontent.com/ale94mleon/BindFlow/main/install/uninstaller.sh)" -- [options]

Removes the BindFlow environments created by the installer. It lists what it will remove and asks
first; it also offers to remove the files Espaloma leaves in your home directory.

  -y, --yes          Do not ask; remove the environments (and the Espaloma files with --caches)
  --caches           Also remove the Espaloma model and the DGL config file
  --dry-run          Only show what would be removed
  --name NAME        Main environment name (default: ${MAIN_ENV})
  --toff-name NAME   Parameterization environment name (default: the one linked to the main
                     environment by the installer, else ${TOFF_ENV})
  -h, --help         Show this help

Requires micromamba or mamba >= 2.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -y|--yes) ASSUME_YES=1 ;;
        --caches) CACHES=1 ;;
        --dry-run) DRY_RUN=1 ;;
        --name) MAIN_ENV="$2"; shift ;;
        --toff-name) TOFF_ENV="$2"; TOFF_ENV_GIVEN=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) usage; printf 'ERROR: Unknown option: %s\n' "$1" >&2; exit 1 ;;
    esac
    shift
done

if [[ -t 1 && -z "${NO_COLOR:-}" ]]; then
    BOLD=$'\033[1m' DIM=$'\033[2m' RED=$'\033[31m' GREEN=$'\033[32m' YELLOW=$'\033[33m' CYAN=$'\033[36m' RESET=$'\033[0m'
else
    BOLD="" DIM="" RED="" GREEN="" YELLOW="" CYAN="" RESET=""
fi

die() {
    printf '%sERROR:%s %s\n' "$RED" "$RESET" "$*" >&2
    exit 1
}
warn() { printf '%sWarning:%s %s\n' "$YELLOW" "$RESET" "$*"; }

# ---------------------------------------------------------------------------
# Conda frontend: the same lookup as the installer ($MAMBA_EXE, micromamba, mamba)
# ---------------------------------------------------------------------------
MM=""
for candidate in "${MAMBA_EXE:-}" "$(command -v micromamba 2>/dev/null || true)" "$(command -v mamba 2>/dev/null || true)"; do
    if [[ -n "$candidate" && -x "$candidate" ]]; then MM="$candidate"; break; fi
done
[[ -n "$MM" ]] || die "micromamba or mamba not found, so the environments cannot be listed or removed."

# Questions go through /dev/tty (as in the installer), so they also work with curl ... | bash.
# Without -y the uninstaller only removes something after a "yes" at the terminal.
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

# env_prefix NAME: full path of environment NAME, empty if it does not exist.
# Read from the environment list, without running anything inside the environment (it may be broken).
env_prefix() {
    "$MM" env list --json | tr ',[]' '\n\n\n' | sed -n 's|^[[:space:]]*"\(.*/envs/'"$1"'\)"[[:space:]]*$|\1|p' | head -n 1
}

# has_package PREFIX NAME: true if pip installed NAME in PREFIX (BindFlow and TOFF come from pip)
has_package() {
    local d
    for d in "$1"/lib/python*/site-packages/"$2"-*.dist-info; do
        if [[ -d "$d" ]]; then return 0; fi
    done
    return 1
}

size_of() { du -sh "$1" 2>/dev/null | cut -f 1; }

printf '\n%sBindFlow uninstaller%s\n\n' "$BOLD" "$RESET"

# ---------------------------------------------------------------------------
# What is installed
# ---------------------------------------------------------------------------
# The installer writes the path of the parameterization environment into the main environment
# (etc/conda/activate.d/bindflow-toff.sh). If it is there, that is the TOFF environment of this
# install, whatever its name.
remove_names=()
remove_prefixes=()

main_prefix="$(env_prefix "$MAIN_ENV")"
link_file="${main_prefix:+$main_prefix/etc/conda/activate.d/bindflow-toff.sh}"
if [[ -z "$TOFF_ENV_GIVEN" && -n "$link_file" && -f "$link_file" ]]; then
    linked="$(sed -n 's/^export BINDFLOW_TOFF_PREFIX="\(.*\)"$/\1/p' "$link_file")"
    if [[ -n "$linked" ]]; then TOFF_ENV="$(basename "$linked")"; fi
fi
toff_prefix="$(env_prefix "$TOFF_ENV")"

# add_env NAME PREFIX PACKAGE: plan to remove the environment if it is really BindFlow's.
# An environment that only shares the name (no BindFlow/TOFF in it) is left alone.
add_env() {
    local name="$1" prefix="$2" package="$3"
    if [[ -z "$prefix" ]]; then
        printf '  %s%s: not found%s\n' "$DIM" "$name" "$RESET"
        return
    fi
    if ! has_package "$prefix" "$package" && [[ ! -f "$prefix/etc/conda/activate.d/bindflow-toff.sh" ]]; then
        warn "'$name' ($prefix) does not contain ${package}; it is not a BindFlow environment and is left alone."
        return
    fi
    if [[ "${CONDA_PREFIX:-}" == "$prefix" ]]; then
        die "'$name' is the active environment. Run 'micromamba deactivate' (or 'mamba deactivate') first."
    fi
    remove_names+=("$name")
    remove_prefixes+=("$prefix")
}
printf '%sEnvironments%s\n' "$BOLD" "$RESET"
add_env "$MAIN_ENV" "$main_prefix" bindflow
add_env "$TOFF_ENV" "$toff_prefix" toff
i=0
while [[ $i -lt ${#remove_names[@]} ]]; do
    printf '  %s  %s%s  (%s)%s\n' "${remove_names[$i]}" "$DIM" "${remove_prefixes[$i]}" "$(size_of "${remove_prefixes[$i]}")" "$RESET"
    i=$((i + 1))
done

# Files outside the environments. Only Espaloma writes any (checked for BindFlow, TOFF and
# openmmforcefields, which write nothing to the home directory):
#   - its model, downloaded by PyTorch into $TORCH_HOME/hub/checkpoints (default ~/.cache/torch).
#     That folder is PyTorch's cache for everything, so only the espaloma-*.pt files are removed.
#   - DGL's config file, written the first time dgl is imported. ~/.dgl can also hold DGL datasets,
#     so only config.json is removed, and the folder only if it is then empty.
cache_files=()
torch_checkpoints="${TORCH_HOME:-${XDG_CACHE_HOME:-$HOME/.cache}/torch}/hub/checkpoints"
for f in "$torch_checkpoints"/espaloma-*.pt; do
    if [[ -f "$f" ]]; then cache_files+=("$f"); fi
done
dgl_dir="${DGLDEFAULTDIR:-$HOME/.dgl}"
if [[ -f "$dgl_dir/config.json" ]]; then cache_files+=("$dgl_dir/config.json"); fi

if [[ ${#cache_files[@]} -gt 0 ]]; then
    printf '\n%sFiles left by Espaloma%s %s(optional; the DGL config is shared with any other DGL install)%s\n' "$BOLD" "$RESET" "$DIM" "$RESET"
    for f in "${cache_files[@]}"; do printf '  %s  %s(%s)%s\n' "$f" "$DIM" "$(size_of "$f")" "$RESET"; done
fi
printf '\n'

if [[ ${#remove_names[@]} -eq 0 && ${#cache_files[@]} -eq 0 ]]; then
    printf 'Nothing to remove.\n'
    exit 0
fi

# ---------------------------------------------------------------------------
# Confirm
# ---------------------------------------------------------------------------
remove_caches=""
if [[ -n "$DRY_RUN" ]]; then
    printf 'Dry run: nothing was removed.\n'
    exit 0
fi
if [[ -n "$INTERACTIVE" ]]; then
    if [[ ${#remove_names[@]} -gt 0 ]]; then
        names="$(printf '%s, ' "${remove_names[@]}")"
        ask_yes_no "Remove ${names%, }?" n "They will be removed" "Kept" || remove_names=()
    fi
    if [[ ${#cache_files[@]} -gt 0 ]]; then
        if [[ -n "$CACHES" ]] || ask_yes_no "Also remove the files left by Espaloma?" n "They will be removed" "Kept"; then
            remove_caches=1
        fi
    fi
elif [[ -n "$ASSUME_YES" ]]; then
    remove_caches="$CACHES"
else
    die "No terminal to ask for confirmation. Run with -y to remove without asking (add --caches for the Espaloma files)."
fi

# ---------------------------------------------------------------------------
# Remove
# ---------------------------------------------------------------------------
# micromamba's output is kept in a temporary file and shown only if a removal fails. Even a clean
# removal prints libmamba warnings and the output of packages' own unlink scripts (e.g. an "rm" of a
# file that does not exist), which look like problems but are not.
removed=""
mm_log="$(mktemp)"
trap 'rm -f "$mm_log"' EXIT
for name in ${remove_names[@]+"${remove_names[@]}"}; do
    printf 'Removing %s ...' "$name"
    if ! "$MM" env remove -n "$name" -y > "$mm_log" 2>&1; then
        printf '\n'
        cat "$mm_log"
        die "Could not remove '$name' (micromamba output above)."
    fi
    printf '\r%s✔%s Removed %s      \n' "$GREEN" "$RESET" "$name"
    removed=1
done
if [[ -n "$remove_caches" ]]; then
    for f in ${cache_files[@]+"${cache_files[@]}"}; do
        rm -f "$f"
        printf '%s✔%s Removed %s\n' "$GREEN" "$RESET" "$f"
    done
    rmdir "$dgl_dir" 2>/dev/null || true   # only if empty
    removed=1
fi
if [[ -z "$removed" ]]; then
    printf 'Nothing was removed.\n'
    exit 0
fi

if [[ ${#remove_names[@]} -eq 0 ]]; then exit 0; fi
printf '\n%sBindFlow was uninstalled.%s\n' "$BOLD" "$RESET"
if [[ -z "$remove_caches" && ${#cache_files[@]} -gt 0 ]]; then
    printf '%sThe Espaloma files listed above were kept; run again with --caches to remove them.%s\n' "$DIM" "$RESET"
fi
# The package cache belongs to micromamba and is shared by all environments, so it is only mentioned
printf '%smicromamba still keeps the downloaded packages in its cache; "micromamba clean --all" frees that space.%s\n' "$DIM" "$RESET"
