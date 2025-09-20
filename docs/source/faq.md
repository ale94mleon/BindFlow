# ❓FAQ

````{dropdown} Atom X in residue Y was not found in rtp entry Y with x atoms while sorting atoms
:color: info
:animate: fade-in-slide-down

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
````

````{dropdown} Intramolecular Interactions: To Couple or Not to Couple
:color: info
:animate: fade-in-slide-down

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

If you use `submit = True` for the functions {py:func}`bindflow.runners.calculate`. One job will be launched to the cluster with the only aim of launching the Snakemake jobs to the cluster and waiting till the completion of the entire workflow. This is inefficient: actively allocated resources in the cluster that are not been used.

A workaround in case a frontend is available is to set `submit = False` and then on the `approach` directory do:

```bash
conda activate BindFlow
cd <approach_directory>
nohup nice -19 ./job.sh > RuleThemAll.out 2>&1 &
```

Now, even if you close your terminal, the process will continue running in the background because of the use of `nohup`. This process is mainly idle, but by using `nice -19`, we lower its priority, so it does not interfere with any main processes running on your front end. You can also use other persistent terminals like [screen](https://www.gnu.org/software/screen/manual/screen.html) or [byobu](https://www.byobu.org/).
````

````{dropdown} Error on fep_ana_get_dg_complex_contributions
:color: info
:animate: fade-in-slide-down

Check {ref}`debugging-bindflow-runs` runs section.

If, in the `.err` file under `slurm_logs`, you find an error such as:


```{error}
Duplicate time values found; it's generally advised to use slicing on DataFrames with unique time values for each row. Use `force=True` to ignore this error.
```

This usually means that BindFlow was restarted and the GROMACS simulation did not handle the restart correctly.

### Step 1 — Identify the problematic window

From the error log, determine which `ligand` and `replica` failed. For example, assume the error points to ligand3 and replica 1

### Step 2 — Detect duplicate time values in XVG files or corrupted XVG files

```python
import numpy as np
from glob import glob 

root_path = "fep/openff_unconstrained-2.0.0/ligand3/1/complex/fep/simulation/"
for xvg in glob(root_path+"/*/prod/prod.xvg"):
    try:
        data = np.loadtxt(xvg, comments=["#", "@"])[:, 0]
    except ValueError:
        print(xvg)
    uniques, counts = np.unique(data, return_counts=True)
    duplicates = uniques[counts > 1]
    if len(duplicates) > 0:
        print(xvg, duplicates)
```

Suppose this identifies:

```bash
fep/openff_unconstrained-2.0.0/ligand3/1/complex/fep/simulation/vdw.0/prod/prod.xvg
```

### Step 3 — Clean the problematic window

Delete all files in that window **except prod.mdp**:

```bash
find fep/openff_unconstrained-2.0.0/ligand3/1/complex/fep/simulation/vdw.0/prod -maxdepth 1 ! -name 'prod.mdp' -type f -delete
```

### Step 4 — Reset and rerun

1. Clean the working directory:

```bash
bindflow clean fep/openff_unconstrained-2.0.0
```

2. Rerun the pipeline.

The affected window (e.g., vdw.0) will be repeated.

### TL;DR

```{mermaid}
flowchart TD
    A[Error detected in slurm_logs .err file] --> B[Identify ligand and replica from error log]
    B --> C[Run Python script to scan XVG files]
    C --> D{Duplicates or corrupted found?}
    D -- Yes --> E[Delete all files except prod.mdp in problematic window]
    E --> F[Run 'bindflow clean']
    F --> G[Rerun pipeline]
    D -- No --> H[Check main BindFlow log for other causes]
```
````
