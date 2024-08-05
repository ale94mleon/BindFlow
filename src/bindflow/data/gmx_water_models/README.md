# Descriptions of files and force fields

All configuration and topology files are sourced from GROMACS force fields, available at [GROMACS GitLab - share/top](https://gitlab.com/gromacs/gromacs/-/tree/main/share/top?ref_type=heads). These files contain topologies and configurations for water models and ions within three force field families: AMBER, CHARMM, and OPLS-AA.

It is assumed that inside the same family, the non-bonded interactions are the same (`epsilon` and `sigma` parameters), which is true for the force fields presented in the GROMACS distribution.

The `ffnonbonded.itp` for each family was taken from:

* AMBER: amber99sb-ildn
* CHARMM: charmm27
* OPLS-AA: oplsaa

These `ffnonbonded.itp` files were modified to retain only the `[ atomtypes ]` section, including only atom types related to the water models and ions. This modification prevents potential conflicts with atom-type definitions from user-provided force fields.

The available force fields and their corresponding configuration files are:

```yaml
amber:
  spc: spc216.gro
  spce: spc216.gro
  tip3p: spc216.gro
  tip4p: tip4p.gro
  tip4pew: tip4p.gro
  tip5p: tip5p.gro
charmm:
  spc: spc216.gro
  spce: spc216.gro
  tip3p: spc216.gro
  tip4p: tip4p.gro
  tip5p: tip5p.gro
  tips3p: spc216.gro
oplsaa:
  spc: spc216.gro
  spce: spc216.gro
  tip3p: spc216.gro
  tip4p: tip4p.gro
  tip4pew: tip4p.gro
  tip5p: tip5p.gro
  tip5pe: tip5p.gro
```
