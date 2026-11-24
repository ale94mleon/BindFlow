# 🗒️ Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.18.0] - 2026.11.24

### Upgrading

- **Slurm:** jobs are now submitted with Snakemake's cluster-generic executor plugin. Existing `BindFlow` environments need it: run the installer again, or `micromamba install -n BindFlow -c conda-forge -c bioconda snakemake-executor-plugin-cluster-generic`. Commands written by hand with `--cluster`/`--cluster-cancel` must change to `--executor cluster-generic --cluster-generic-submit-cmd 'sbatch --parsable ...' --cluster-generic-cancel-cmd scancel` (see the deploy guide).
- **Environment modules:** if `module: command not found` appears on compute nodes, load the module system as the first entry of `extra_directives: dependencies`, e.g. `source /etc/profile.d/modules.sh` (or `lmod.sh`, depending on the cluster). BindFlow runs those commands in a non-login shell, where `module` is only defined if it was exported.
- **Python:** 3.11 or newer is required (Snakemake 9 needs it). The installer already uses Python 3.11.
- **pip:** `pip install bindflow` installs BindFlow alone, without dependencies. Install with the installer; inside an existing environment, `pip install --upgrade bindflow` updates BindFlow without touching the conda packages.
- **`bindflow clean`** no longer kills every Snakemake process or cancels all of your Slurm jobs; it only acts on the given pipeline folder. Use `--dry-run` to see what it would do.

### Added

- The installer runs straight from GitHub (`/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/ale94mleon/BindFlow/main/install/installer.sh)"`). In a terminal, it asks whether to add MM(PB/GB)SA and Espaloma and whether to recreate existing environments, and confirms before installing; `-y` skips the questions (`--help` lists all options). The installation docs now show this single command.
- `install/uninstaller.sh`: removes the `BindFlow` and `BindFlow-TOFF` environments after listing them and asking, and optionally the files Espaloma leaves in the home directory (its model and DGL's config). It only removes environments that contain BindFlow or TOFF; `--dry-run`, `-y` and `--caches` are available. See the installation docs.
- The installer stops on native Windows (Git Bash, MSYS2, Cygwin) and points to WSL2, where BindFlow installs as on Linux. AmberTools and GROMACS have no Windows builds on conda-forge.
- `bindflow clean --dry-run` (`-n`): lists the processes, Slurm jobs and folders that would be affected, without touching them.
- When it finishes, the installer lists the next steps: install GROMACS (with a link to the docs) and run `test_fep` (and `test_mmpbsa` with `--mmpbsa`) to check that BindFlow runs simulations with it.

### Changed

- The installer takes the environment files from the local clone only when it runs as a file or with `--editable`; otherwise it downloads them (before, it looked in the current directory).
- Slurm jobs are submitted through the [cluster-generic executor plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/cluster-generic.html) (`--executor cluster-generic`), now part of the `BindFlow` environment (see Upgrading). The deploy guide and the A2A example show the new command.
- `bindflow clean` only acts on the given pipeline: it cancels the Slurm jobs whose working directory is inside the pipeline folder and stops the Snakemake processes running there (and the jobs they started). Before, it killed every Snakemake process and cancelled all of the user's Slurm jobs.
- The installer and uninstaller only run once fully downloaded: a transfer cut off midway no longer runs part of the script.
- Python >= 3.11 is required: Snakemake 9 and its cluster-generic plugin need it, and the code did not run on 3.8 and 3.9.
- BindFlow no longer declares dependencies for pip (`pyproject.toml` had an incomplete list). They come from the environment files in `install/envs`, and several are not on PyPI. `pip install bindflow` installs BindFlow alone, so it can no longer replace conda packages when `--no-deps` is forgotten. If a dependency is missing, BindFlow now says which one and points to the installation docs, at import time or before submitting anything (the `snakemake` command, the cluster-generic plugin for Slurm, MDRestraintsGenerator for FEP, gmx_MMPBSA for MM(PB/GB)SA).
- {py:func}`bindflow.runners.calculate` raises a `ValueError` when two ligand files have the same name (e.g. `a/lig1.sdf` and `b/lig1.sdf`): the name is the ligand's output folder, so their results overwrote each other.

### Fixed

- Slurm runs failed at once with `unrecognized arguments`: the generated `job.sh` used `--cluster` and `--cluster-cancel`, which Snakemake 8 removed (BindFlow pins Snakemake 9.24).
- Boolean Slurm options (e.g. `exclusive: true`) were passed as `--exclusive=`, which `sbatch` rejects; they are now passed as `--exclusive`.
- The gmx_MMPBSA step always failed when no `dependencies` were configured: the command started with a bare `&&`.
- FEP analysis asked for `prod.mdp` files of windows that do not exist when the number of vdw/coul lambda windows differed from the complex bonded ones.
- `bindflow clean` crashed on subdirectories inside `slurm_logs`.
- {py:func}`bindflow.utils.tools.gmx_runner` raised `TypeError` with the default `load_dependencies=None`.
- {py:func}`bindflow.utils.tools.unarchive` could write outside the target folder when given a crafted archive (`../` paths, links). Extraction now uses tarfile's `data` filter (new {py:func}`bindflow.utils.tools.safe_extract` and {py:func}`bindflow.utils.tools.safe_extract_all`).
- {py:func}`bindflow.utils.tools.archive` stored XTC files under their absolute path, so {py:func}`bindflow.utils.tools.unarchive` restored them under `target/home/...` instead of inside the project.
- mdrun options set to `0` (e.g. `gpu_id: 0`) were silently dropped.
- A custom `Slipids_2020` force field (`custom_ff_path`) was always shadowed by the bundled copy.
- Converting water `[ settles ]` to `[ constraints ]` wrote the topology twice when the molecule had no `[ settles ]` section; it now raises a clear error.
- `MDP(**parameters)` changed the default parameters of every MDP created afterwards.
- {py:func}`bindflow.free_energy.gather_results.get_raw_fep_data` cut ligand names at the last dot (`lig.1` became `lig`).

## [0.17.0] - 2026.11.21

### Upgrading

Small molecules are now parameterized in a separate `BindFlow-TOFF` environment. Updating BindFlow with `pip` in an existing environment is not enough: BindFlow stops before submitting anything if it cannot find `BindFlow-TOFF`. Run the installer (`install/installer.sh`) to create both environments, or pass `toff_prefix` to {py:func}`bindflow.runners.calculate`. The Python API and the configuration files are otherwise unchanged.

### Added

- `install/installer.sh`: installs BindFlow in two linked environments created with micromamba (or mamba >= 2): `BindFlow` (workflow, analysis and, with `--mmpbsa`, MM(PB/GB)SA) and `BindFlow-TOFF` (small-molecule parameterization and pdbfixer; Espaloma with `--espaloma`). The environment files are in `install/envs`.
- `bindflow.utils.toff_env`: runs TOFF and pdbfixer from the `BindFlow-TOFF` environment, in a separate process.
- Keyword `toff_prefix` to {py:func}`bindflow.runners.calculate`. If not given, `BINDFLOW_TOFF_PREFIX` (set by the installer) or an environment called `BindFlow-TOFF` next to the current one is used. Before submitting anything, BindFlow checks that the environment exists and has TOFF >= 0.3.0, when it is needed.
- Development container (`.devcontainer`) built with the installer, including GROMACS 2025.4 from conda-forge for the tests.
- GitHub Actions workflow (`test-installation.yml`) that tests the installer on Linux (x86_64, ARM) and macOS Apple Silicon, on changes to the installer or the code. It runs the small-molecule tests only; `test_fep` and `test_mmpbsa` are too expensive for CI and are run locally.

### Changed

- TOFF, OpenMM, OpenFF, RDKit and pdbfixer are no longer needed in the main BindFlow environment.
- Python 3.11 is supported.
- GROMACS 2026 is supported (range [2022, 2027)).
- Installation docs rewritten around the installer.

### Removed

- Support for macOS Intel (x86_64): the installer stops there. conda-forge's `antechamber` on osx-64 randomly reads zero atoms from SDF files, so OpenFF and GAFF parameterization fail intermittently (upstream bug, see TOFF's known issues).
- The interactive `environment.yml` builder of the docs (`custom.js`, `custom.css`) and the old environment files in `docs/source/env`.

## [0.16.3] - 2026.05.10

### Fixed

- [#52](https://github.com/ale94mleon/BindFlow/issues/52): Mismatch between tests and package.
- Update docs.
- `preparation.system_builder.system_combiner`. Use slicing instead deep copy so the Structure class initialize properly avoiding errors such as `AttributeError: 'GromacsTopologyFile' object has no attribute 'symmetry'`.

### Added

- Keyword `cofactor_selection`, this is a GMX selection for better control of the definition of complex cofactors given via .gro/.top files. For example: "resname GDP or resname GTP or resname MG". This is used for the definition of the groups in the thermostat for membrane systems, This is only relevant for membrane systems. For soluble systems everything is treated as a single group.
- `solv_d`, `solv_bt` (`d` and `bt` flags of `gmx editconf`), `solv_rmin`, and `solv_ion_conc` (`rmin` and `conc` flags of `gmx genion`) to {py:func}`bindflow.runners.calculate`. Now there is more room to control the solvation step.

### Changed

- Tests are part of the wheel, `pytest --pyargs bindflow.tests` is used for testing.
- The box type for soluble complex and the ligand simulation was hard coded to `octahedron`, now it is customizable with `solv_bt`, but the default value is the more efficient `dodecahedron` box.

## [0.15.1] - 2026.03.17

### Added

- CI/CD for PyPI deployment
- Check for `gmx` version

### Changed

- Update documentation, tutorials and installation
- Improve code styling
- Use always 1 bar instead 1.01325 for membrane simulations during equilibration

### Fixed

- Better logging
- Move to new setuptools and update versioning
- Exclude `pme` flag during equilibration
- Pass `gmx` dependencies during solvation

### Removed

- MDP options from `bonded` and `coul` windows as there are not being used. Those flags are only relevant when either Lennard-Jones or Coulomb interaction are decoupled. In our case for `coul` they also irrelevant as the options are applied only if `sc-coul = yes`, which by default is `no`. Check the [GROMACS docs](https://manual.gromacs.org/current/user-guide/mdp-options.html#mdp-sc-coul) for more info on these MDP options.

```mdp
sc-alpha                 = 0.5
sc-power                 = 1
sc-sigma                 = 0.3
```

- `"MDRestraintsGenerator @ git+https://github.com/ale94mleon/MDRestraintsGenerator.git@dev"` dependency from `pyproject.toml` as it is not accepted by PyPi. This dependency is built during the conda environment step anyway.

[unreleased]: https://github.com/ale94mleon/bindflow/compare/v0.17.0...HEAD
[0.17.0]: https://github.com/ale94mleon/bindflow/compare/v0.16.3...v0.17.0
[0.16.3]: https://github.com/ale94mleon/bindflow/compare/v0.15.1...v0.16.3
[0.15.1]: https://github.com/ale94mleon/bindflow/compare/v0.13.0...v0.15.1
