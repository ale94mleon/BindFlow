# FAQ

<!-- :::{dropdown} Atom X in residue XXX was not found in rtp entry XX with X atoms while sorting atoms
:open: -->

## Atom X in residue XXX was not found in rtp entry XX with X atoms while sorting atoms

BindFlow uses [PDBFixer](https://github.com/openmm/pdbfixer) for the standardization of PDB files. However, [PDBFixer](https://github.com/openmm/pdbfixer) is not bulletproof. For such cases where [pdb2gmx](https://manual.gromacs.org/current/onlinehelp/gmx-pdb2gmx.html) fails with the generated PDB file from PDBFixer; the user may fix the PDB by hand.

Another strategy is give to BindFlow the TOP and GRO files instead of the plain PDB file.

### Example

```bash
wget https://raw.githubusercontent.com/openforcefield/protein-ligand-benchmark/main/data/mcl1/01_protein/crd/protein.pdb
```

The following code will fail

```python
from bindflow.preparation.system_builder import MakeInputs

ligands = [
    "lig_63.mol",
]

protein = {
    'conf': 'protein.pdb',
}

builder = MakeInputs(
    protein=protein,
    water_model='amber/tip3p',
    hmr_factor=2.5,
    builder_dir="builder")
builder(ligand_definition=ligands[0], out_dir='test_test')
```

The error is going to be:

```txt
-------------------------------------------------------
Program:     gmx pdb2gmx, version 2023.3-dev-20231019-5e5ea27-local
Source file: src/gromacs/gmxpreprocess/pdb2gmx.cpp (line 870)

Fatal error:
Atom C in residue NME 322 was not found in rtp entry NME with 6 atoms
while sorting atoms.
.

For more information and tips for troubleshooting, please check the GROMACS
website at http://www.gromacs.org/Documentation/Errors
-------------------------------------------------------
```

#### Getting TOP and GRO files and changing protein definition

```python
import parmed
import os
from pathlib import Path
from openmm import app

pdb_prot_file = Path("inputs/protein/protein.pdb")


# Any force field in
# https://github.com/openmm/openmm/tree/master/wrappers/python/openmm/app/data
# can be selected

protein_force_field = 'amber14-all'
forcefield = app.ForceField(f'{protein_force_field}.xml', 'tip3p.xml')

dirname = pdb_prot_file.parent
filename = pdb_prot_file.stem
out_dir = dirname / f"protein-{protein_force_field}"
Path(out_dir).mkdir(exist_ok=True)
print(out_dir)

pdb_obj = app.PDBFile(str(pdb_prot_file))
openmm_topology = pdb_obj.topology
# Create an OpenMM System from an OpenMM Topology object
system = forcefield.createSystem(pdb_obj.topology)

struct = parmed.openmm.load_topology(pdb_obj.topology, system=system, xyz=pdb_obj.positions)
for file_type in [out_dir / f'{filename}.gro', out_dir / f'{filename}.top']:
    struct.save(str(file_type), overwrite=True)
```

And then use the GRO and TOP files:

```python
protein = {
    'conf': 'inputs/protein/protein-amber14-all/protein.gro',
    'top': 'inputs/protein/protein-amber14-all/protein.top',
}
```

## Intramolecular Interactions: To Couple or Not to Couple

The thermodynamic cycle employed by BindFlow involves the parameter `couple-intramol = yes`, indicating that the intramolecular interactions of the ligand change alongside the λ parameter. For instance, upon removing Coulomb and van der Waals interactions in the water box, the ligand ceases to interact with the solvent while also disengaging from self-interactions via non-bonded interactions.

This configuration proves advantageous for larger molecules wherein intramolecular interactions may occur over considerable distances. Otherwise, distant regions of the molecule would overly interact via explicit pair interactions, leading to artificially strong bonding, which could bias the resulting free energy.

Although `couple-intramol = yes` is primarily beneficial for large molecules, BindFlow sets the default value to `yes` due to the uncertainty of the molecules being processed.

Consequently, the energy contributions term `{complex, ligand}_{coul, vdw}` no longer solely denotes the change in free energy for the coupling/decoupling of `{coul, vdw}` interactions of the ligand in the `{protein, solvent}` environment. It encompasses the free energy difference for coupling/decoupling these interactions in the `{protein, solvent}` environment plus the free energy variation when the ligand's intramolecular interaction (`{coul, vdw}`) is also coupled/decoupled. This explains some big values that are usually obtained for the `coul` contribution whcih is mainly formed by the intramolecular contribution

BindFlow uses the same λ-schedule for the simulations of the ligand in the water solvent and in the binding pocket of the protein (same simulation time and same λ-values). This means that the contribution of coupling/decoupling the ligand intramolecular interaction is the same but with the opposite sign for each contribution (either `coul` or `vdw`). This helps us to recover some useful information from the energetic contributions. the following sums will cancel out the ligand intramolecular contribution

- {math}`\text{complex_vdw} + \text{ligand_vdw} = \text{vdw_contrib}` -> estimation of the van der Waal contribution to the binding
- {math}`\text{complex_coul} + \text{ligand_coul} = \text{coul_contrib}` -> estimation of the Coulomb contribution to the binding

Then the equation of the cycle:

```{math}
\Delta G = \text{ligand_coul} + \text{ligand_vdw} - \text{boresch} + \text{complex_vdw} + \text{complex_coul} - \text{bonded}
```

We can rewrite it as:

```{math}
\Delta G = \text{vdw_contrib} + \text{coul_contrib} - \text{boresch} - \text{bonded}
```

## Where to run the main job?

If you use `submit = True` for the functions `bindflow.calculate_abfe` or `bindflow.calculate_calculate_mmpbsa`. One job will be launched to the cluster with the only aim of launching the Snakemake jobs to the cluster and waiting till the completion of the entire workflow. This is inefficient: actively allocated resources in the cluster that are not been used.

A workaround in case a frontend is available is to set `submit = False` and then on the `approach` directory do:


```bash
conda activate BindFlow
```

## Cancelling running jobs in a Slurm HPC cluster