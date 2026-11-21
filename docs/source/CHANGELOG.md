# 🗒️ Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]

### Added

- The installer runs straight from GitHub (`/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/ale94mleon/BindFlow/main/install/installer.sh)"`). In a terminal, it asks whether to add MM(PB/GB)SA and Espaloma and whether to recreate existing environments, and confirms before installing; `-y` skips the questions (`--help` lists all options). The installation docs now show this single command.
- The installer stops on native Windows (Git Bash, MSYS2, Cygwin) and points to WSL2, where BindFlow installs as on Linux. AmberTools and GROMACS have no Windows builds on conda-forge.

### Changed

- The installer takes the environment files from the local clone only when it runs as a file or with `--editable`; otherwise it downloads them (before, it looked in the current directory).

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
