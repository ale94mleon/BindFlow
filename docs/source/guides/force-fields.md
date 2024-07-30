# Force fields

BindFlow provides several out-of-the-box force fields. This section explains how to access them and integrate custom force fields within BindFlow. We mainly focus on the input possibilities for [BindFlow's runner](#bindflow-runners) functions.

BindFlow offers a variety of force field options, but as Uncle Ben says, "With great power comes great responsibility." Users must choose the appropriate combination of force fields. By default, BindFlow offers a suitable combination.

## Structure inputs

Six keywords control the type of force field used for each specific component in the system:

1. `protein`: Definition of the host
2. `membrane`: Definition of the membrane
3. `ligands`: A list of ligand's definitions
4. `cofactor`: Definition of the cofactor
5. `water_model`: Type of GROMACS' water model to use
6. `custom_ff_path`: Path to the custom force field if needed.

## Partial and Full Definitions

For a straightforward setup, you can provide the path to the corresponding file(s), which we will call the _partial definition_. However, you also have the option to fine-tune the definition of force fields for each component, referred to as the _full definition_.

In the following examples, we will use the runner {py:func}`bindflow.run_abfe.calculate_abfe`, but the same applies to {py:func}`bindflow.run_mmpbsa.calculate_mmpbsa`.

````````{tab} Partial definition
``````{tab} protein
The protein will be processed with [amber99sb-ildn](https://ambermd.org/#ff) force field after been fixed with [PDBFixer](https://github.com/openmm/pdbfixer).

```{hint}
It is advised to spend some time on the processing of the protein beforehand (better a minute than repeat the whole campaign):

1. Missing atoms
2. Missing loops
3. Terminal capping
4. Protonation state

All the above steps are highly system-dependent, and while PDBFixer can handle some minor issues, it is far from perfect. In addition, our use of PDBFixer is very simple.
```

```python
calculate_abfe(
    ...
    protein="path/to/protein.{pdb;gro}",
    ...   
)
```
``````
``````{tab} membrane

The membrane will be processed with [SLipids_2020](http://www.fos.su.se/~sasha/SLipids/Downloads.html).

```{dropdown} Getting the membrane.pdb file
:color: info
:animate: fade-in-slide-down
:icon: rocket

For a membrane systems you first need to embed the protein into the membrane. This can easily done with [CHARMM-GUI](https://www.charmm-gui.org). This is a non-exhaustive list of steps for this process:

Processing on [CHARMM-GUI](https://www.charmm-gui.org):

* ACE and CT3 terminus
* pH=7
* Run PPM 2.0
* Hexagonal box
* Water thickness: 15
* Length of X and Y: 90 (initial guess)
* Only POPC for a simple membrane
* Do not include ions
* force field:
    - AMBER:
        - Protein: ff14sb
        - Lipid: SLipids
        - Water: TIP3P
* Input Generation Options: GROMACS


Open `gromac/step5_input.gro`. This file has the crystal information which is important during the solvation step in BindFLow. Convert to PDB with `gmx editconf -f step5_input.gro -o temporal.pdb `. The last PDB will also have the crystal info. Split the PDB in POPC and protein.

In [PyMOL](https://www.pymol.org):

* `select popc, resn POPC`
* `select prot, (polymer.protein or resn ACE or resn NME)`

👀 Depending on how you are processing the protein, you might need to change the entry ATOM to HETATM for the CAP groups ACE and NME manually in the PDB file.
```

```python
calculate_abfe(
    ...
    membrane="path/to/membrane.pdb",
    ...   
)
```
``````
``````{tab} ligands
```python
calculate_abfe(
    ...
    ligands=[
        "path/to/ligand1.{mol;sdf}",
        "path/to/ligand2.{mol;sdf}",
        "path/to/ligand3.{mol;sdf}",
        ...
    ],
    ...   
)
```
``````
``````{tab} cofactor
```python
calculate_abfe(
    ...
    cofactor="path/to/cofactor.{mol;sdf}",
    ...   
)
```
``````
````````

````````{tab} Full definition
``````{tab} protein
`````{tab} by code
````{tab} on GROMACS distribution

You can access all the [GROMACS force fields](https://manual.gromacs.org/current/user-guide/force-fields.html) by their code, they will be pass to [pdb2gmx](https://manual.gromacs.org/documentation/current/onlinehelp/gmx-pdb2gmx.html) through the flag `-ff` after been fixed with [PDBFixer](https://github.com/openmm/pdbfixer).

```{hint}
It is advised to spend some time on the processing of the protein beforehand (better a minute than repeat the whole campaign):

1. Missing atoms
2. Missing loops
3. Terminal capping
4. Protonation state

All the above steps are highly system-dependent, and while PDBFixer can handle some minor issues, it is far from perfect. In addition, our use of PDBFixer is very simple
```

```python
calculate_abfe(
    ...
    protein={
        "conf": "path/to/protein.{pdb;gro}",
        "ff":{
            "code": <GMX_force_field_code>
        }
    },
    ...   
)
```
````
````{tab} external

To add even more flexibility, you can use any external force field ported to GROMACS, in this case you just need to copy your `force_field.ff` (e.g. `charmm36-jul2022.ff`) to your desired directory and pass the path of this directory to `custom_ff_path` parameter. If you have more force fields, you can copy all of them in the same directory. BindFlow will internally set the following environmental variable at run time.

```python
os.environ["GMXLIB"] = os.path.abspath(custom_ff_path)
```

```{warning}
See how the force field directory ends in `.ff`; e.g. `charmm36-jul2022.ff`. This is needed.
```

The force field code (e.g. for `charmm36-jul2022.ff`, the code is `charmm36-jul2022`) will be pass to [pdb2gmx](https://manual.gromacs.org/documentation/current/onlinehelp/gmx-pdb2gmx.html) through the flag `-ff` after been fixed with [PDBFixer](https://github.com/openmm/pdbfixer).

```{hint}
It is advised to spend some time on the processing of the protein beforehand (better a minute than repeat the whole campaign):

1. Missing atoms
2. Missing loops
3. Terminal capping
4. Protonation state

All the above steps are highly system-dependent, and while PDBFixer can handle some minor issues, it is far from perfect. In addition, our use of PDBFixer is very simple
```

```python
calculate_abfe(
    ...
    protein={
        "conf": "path/to/protein.{pdb;gro}",
        "ff":{
            "code": <external_force_field_code>
        }
    },
    custom_ff_path='parent/directory/of/custom.ff'
    ...   
)
```
````
`````
`````{tab} by top

```{admonition} Be careful
:class: danger

It is advised to build a single topology file without any `include` statements. If you want to use those include statements, they **MUST BE** absolute paths to their corresponded files.
```

```python
calculate_abfe(
    ...
    protein={
        "conf": "path/to/protein.gro",
        "top": "path/to/protein.top",
    },
    ...   
)
```
`````
``````
``````{tab} membrane

In this case a topology must be generated and provided. This topology can also be obtained from CHARMM-GUI.

```{admonition} Be careful
:class: danger

It is advised to build a single topology file without any `include` statements. If you want to use those include statements, they **MUST BE** absolute paths to their corresponded files.

You must always past the PDB despite passing a topology, GRO files are not accepted at the moment.
```

```{dropdown} Getting the membrane.pdb file
:color: info
:animate: fade-in-slide-down
:icon: rocket

For a membrane systems you first need to embed the protein into the membrane. This can easily done with [CHARMM-GUI](https://www.charmm-gui.org). This is a non-exhaustive list of steps for this process:

Processing on [CHARMM-GUI](https://www.charmm-gui.org):

* ACE and CT3 terminus
* pH=7
* Run PPM 2.0
* Hexagonal box
* Water thickness: 15
* Length of X and Y: 90 (initial guess)
* Only POPC for a simple membrane
* Do not include ions
* force field:
    - AMBER:
        - Protein: <protein force field>
        - Lipid: <membrane force field>
        - Water: TIP3P
* Input Generation Options: GROMACS


Open `gromac/step5_input.gro`. This file has the crystal information which is important during the solvation step in BindFLow. Convert to PDB with `gmx editconf -f step5_input.gro -o temporal.pdb `. The last PDB will also have the crystal info. Split the PDB in POPC and protein.

In [PyMOL](https://www.pymol.org):

* `select popc, resn POPC`
* `select prot, (polymer.protein or resn ACE or resn NME)`

👀 Depending on how you are processing the protein, you might need to change the entry ATOM to HETATM for the CAP groups.

You may also need to manually split the topology into protein and membrane.
```

```python
calculate_abfe(
    ...
    membrane={
        "conf": "path/to/membrane.pdb",
        "top": "path/to/membrane.top",
    }
    ...
)
```
``````
``````{tab} ligands
`````{tab} openff

Any force field from [OpenFF](https://openforcefield.org/force-fields/force-fields/) can be acceded by setting its name as `code` (see that the `.offxml` extension is kept).

If `code` is not provided, the default force field for `type = "openff"` is `openff_unconstrained-2.0.0.offxml`.

```python
ligand_files = [
    "path/to/ligand1.{mol;sdf}",
    "path/to/ligand2.{mol;sdf}",
    "path/to/ligand3.{mol;sdf}",
    ...
]

ligands = []
for ligand_file in ligand_files:
    ligands.append({
        "conf": ligand_file,
        "ff":{
            "type": "openff",
            "code": "openff_unconstrained-2.0.0.offxml"
        }
    })

calculate_abfe(
    ...
    ligands=ligands,
    ...
)
```
`````
`````{tab} espaloma

It is recommended to use `espaloma >= 0.3.1`.

If `code` is not provided, the default force field for `type = "espaloma"` is `espaloma-0.3.1`.

```python
ligand_files = [
    "path/to/ligand1.{mol;sdf}",
    "path/to/ligand2.{mol;sdf}",
    "path/to/ligand3.{mol;sdf}",
    ...
]

ligands = []
for ligand_file in ligand_files:
    ligands.append({
        "conf": ligand_file,
        "ff":{
            "type": "espaloma",
            "code": "espaloma-0.3.1"
        }
    })

calculate_abfe(
    ...
    ligands=ligands,
    ...
)
```
`````
`````{tab} gaff

If `code` is not provided, the default force field for `type = "gaff"` is `gaff-2.11`.

```python
ligand_files = [
    "path/to/ligand1.{mol;sdf}",
    "path/to/ligand2.{mol;sdf}",
    "path/to/ligand3.{mol;sdf}",
    ...
]

ligands = []
for ligand_file in ligand_files:
    ligands.append({
        "conf": ligand_file,
        "ff":{
            "type": "gaff",
            "code": "gaff-2.11"
        }
    })

calculate_abfe(
    ...
    ligands=ligands,
    ...
)
```
`````
`````{tab} custom force field

```{admonition} Be careful
:class: danger

It is advised to build a single topology file without any `include` statements. If you want to use those include statements, they **MUST BE** absolute paths to their corresponded files.
```

```python
ligands = [
    {
        "conf": "path/to/ligand1.gro",
        "top": "path/to/ligand1.top"
    },
    {
        "conf": "path/to/ligand2.gro",
        "top": "path/to/ligand2.top"
    },
    {
        "conf": "path/to/ligand3.gro",
        "top": "path/to/ligand3.top"
    },
    ...
]

calculate_abfe(
    ...
    ligands=ligands,    
    ...
)
```
`````
``````
``````{tab} cofactor
`````{tab} openff

Any force field from [OpenFF](https://openforcefield.org/force-fields/force-fields/) can be acceded by setting its name as `code` (see that the `.offxml` extension is kept).

If `code` is not provided, the default force field for `type = "openff"` is `openff_unconstrained-2.0.0.offxml`.

```python
calculate_abfe(
    ...
    cofactor={
        "conf": "path/to/cofactor.{mol;sdf}",
        "ff":{
            "type": "openff",
            "code": "openff_unconstrained-2.0.0.offxml"
        }
    },
    ...
)
```
`````
`````{tab} espaloma

It is recommended to use `espaloma >= 0.3.1`.

If `code` is not provided, the default force field for `type = "espaloma"` is `espaloma-0.3.1`.

```python
calculate_abfe(
    ...
    cofactor={
        "conf": "path/to/cofactor.{mol;sdf}",
        "ff":{
            "type": "espaloma",
            "code": "espaloma-0.3.1"
        }
    },
    ...
)
```
`````
`````{tab} gaff

If `code` is not provided, the default force field for `type = "gaff"` is `gaff-2.11`.

```python
calculate_abfe(
    ...
    cofactor={
        "conf": "path/to/cofactor.{mol;sdf}",
        "ff":{
            "type": "gaff",
            "code": "gaff-2.11"
        }
    },
    ...
)
```
`````
`````{tab} custom force field
````{tab} non-water

```{admonition} Be careful
:class: danger

It is advised to build a single topology file without any `include` statements. If you want to use those include statements, they **MUST BE** absolute paths to their corresponded files.
```

```python
calculate_abfe(
    ...
    cofactor={
        "conf": "path/to/cofactor.gro",
        "top": "path/to/cofactor.top",
    },    
    ...
)
```
````
````{tab} water-like

```{admonition} Be careful
:class: danger

It is advised to build a single topology file without any `include` statements. If you want to use those include statements, they **MUST BE** absolute paths to their corresponded files.

In the case that the cofactor(s) is (are) water-like molecule(s), this should be specified by the keyword `is_water = True`. In this case, a special treatment is done in BindFlow internally. Here, its settles section (if any) will be changed to TIP3P-like triangular constraints. Check the discussion [How to treat specific water molecules as ligand?](https://gromacs.bioexcel.eu/t/how-to-treat-specific-water-molecules-as-ligand/3470/9). Note that this is only possible for TIP3P-like water molecules at the moment.
```

```python
calculate_abfe(
    ...
    cofactor={
        "conf": "path/to/cofactor.gro",
        "top": "path/to/cofactor.top",
        "is_water": True,
    },    
    ...
)
```
````
`````
``````
````````

## Water models

BindFlow comes with (at present 07.2024) all water models distributed with GROMACS. They are set by the keyword: `water_model`. E.g.:

```python
calculate_abfe(
    ...
    water_model="amber/tip3p"
    ...
)
```

The structure of the string is `force_field_family/water_model`; `amber/tip3p` is the default. See {py:meth}`bindflow.preparation.solvent.Solvate.__init__` for a list of all available water models.

## MDP options modification based on the force field. AMBER and CHARMM36-like force fields example

Some Molecular Dynamic Parameters (MDP) are usually, rather than interchangeable options, parts of each force field derivation and parametrization. So, we should always use those parameters during our simulations. A typical example are AMBER and CHARMM36-like force fields:

````{tab} AMBER-like force fields

BindFlow use this parameter by default. So, you do not need to modify them.

```yaml
constraints: all-bonds
cutoff-scheme: Verlet
vdwtype: cutoff
vdw-modifier: Potential-Shift-Verlet
rlist: 1.2
rvdw: 1.0
coulombtype: PME
rcoulomb: 1.0
DispCorr: EnerPres
```
````

````{tab} CHARMM36-like force fields

In this case you should pass to to all the steps these parameters. BindFlow works with AMBER-like force fields by default

```yaml
constraints: h-bonds
cutoff-scheme: Verlet
vdwtype: cutoff
vdw-modifier: force-switch
rlist: 1.2
rvdw: 1.2
rvdw-switch: 1.0
coulombtype: PME
rcoulomb: 1.2
DispCorr: no
```
````

## Final note

Remember to cite properly the main references if you use any of the force fields in your work.
