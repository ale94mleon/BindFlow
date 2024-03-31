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

```
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
os.makedirs(out_dir, exist_ok=True)
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
<!-- ::: -->
