# 🗒️ Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]

### Fixed

- Update docs.
- `preparation.system_builder.system_combiner`. Use slicing instead deep copy so the Structure class initialize properly avoiding errors such as `AttributeError: 'GromacsTopologyFile' object has no attribute 'symmetry'`.

### Added

- Keyword `cofactor_selection`, this is a GMX selection for better control of the definition of complex cofactors given via .gro/.top files. For example: "resname GDP or resname GTP or resname MG". This is used for the definition of the groups in the thermostat for membrane systems, This is only relevant for membrane systems. For soluble systems everything is treated as a single group.
- `solv_d`, `solv_bt` (`d` and `bt` flags of `gmx editconf`), `solv_rmin`, and `solv_ion_conc` (`rmin` and `conc` flags of `gmx genion`) to {py:func}`bindflow.runners.calculate`. Now there is more room to control the solvation step.

### Changed

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

[unreleased]: https://github.com/ale94mleon/bindflow/compare/v0.15.1...HEAD
[0.15.1]: https://github.com/ale94mleon/bindflow/compare/v0.13.0...v0.15.1
