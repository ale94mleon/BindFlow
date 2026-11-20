#!/usr/bin/env bash
# Dev environment: BindFlow + BindFlow-TOFF from install/installer.sh (editable, with MM/PBSA),
# plus what is only needed for development: GROMACS (tests), pytest and the docs toolchain.
set -euo pipefail
cd "$(dirname "$0")/.."

# safe.directory avoids git "dubious ownership" errors, which versioningit would hit on the mounted repo.
git config --global --add safe.directory '*'

# Espaloma needs dgl, which conda-forge does not build for linux-aarch64 (e.g. Docker on Apple Silicon)
extra=()
[[ "$(uname -m)" == "x86_64" ]] && extra+=(--espaloma)

bash install/installer.sh --editable . --mmpbsa --force ${extra[@]+"${extra[@]}"}

# Usually building GROMACS for your architecture translate into better computer efficient
# Here we are using the conda version jsut for convinience for the test suit.
micromamba install -n BindFlow -c conda-forge --override-channels -y "gromacs=2025.4" pytest
micromamba run -n BindFlow python -m pip install -r requirements_docs.txt
