# Customizing the Workflow

BindFlow is highly customizable. In this section, we will discuss the options accessible through the `global_config` keyword of the [BindFlow's runners](#bindflow-runners). The `global_config` is a nested Python dictionary. A useful tip is to write this dictionary in YAML format, which is essentially a "human-readable dictionary." Here, we will go through each section of this YAML file:

````{dropdown} config.yml
:color: info
:animate: fade-in-slide-down
:icon: rocket
```yaml
cluster:
  options: # Depending on the Scheduler
    calculation:
      <cluster_options_for_calculation_jobs>
    job: # Optional
      <cluster_options_for_launcher_job>
extra_directives: # Optional
  dependencies:
    - <dependency_commands>
  mdrun:
    ligand:
      <mdrun_keywords_for_ligand_simulation>
    complex:
      <mdrun_keywords_for_complex_simulation>
    all:
      <mdrun_keywords_for_ligand_and_complex_simulation>      
samples: <number_of_samples_for_mmpbsa> # Optional
nwindows: # Optional
  ligand:
      vdw: <number_of_vdw_windows>
      coul: <number_of_coul_windows>
  complex:
      vdw: <number_of_vdw_windows>
      coul: <number_of_coul_windows>
      bonded: <<number_of_bonded_windows>
mmpbsa: # Optional
  <mm(pb/gb)sa_options>
mdp:  # Optional
  ligand:
    equi:
      <step>:
        <mdp_ligand_equi_step_options>
    fep:
      vdw:
        <step>:
          <mdp_ligand_fep_vdw_step_options>
      coul:
        <step>:
          <mdp_ligand_fep_coul_step_options>
  complex:
    equi:
      <step>:
        <mdp_complex_equi_step_options>
    fep:
      vdw:
        <step>:
          <mdp_complex_fep_vdw_step_options>
      coul:
        <step>:
          <mdp_complex_fep_coul_step_options>
      bonded:
        <step>:
          <mdp_complex_fep_bonded_step_options>
    mmpbsa:
      prod:
        <mdp_complex_mm(pb/gb)sa_options>
```
````

```{important}

The options provided on this configuration have higher priority to those ones passed as keywords to the runner function. For example:

* Pass to the runner `dt_max=0.003` and specify for an specific simulation `dt=0.004` (see [MDP section](#mdp-optional))
* Pass to the runner `threads=12` and specify for `mdrun` the flag `-nt 10` (see [mdrun section](#mdrun-optional))

```

## `cluster`

This section specifies the computational resources required by the selected scheduler (see the [BindFlow’s deploy](#bindflow-deploy) section). It allows you to control the allocated resources for two types of jobs:

1. **Main Job (`job`):** This job primarily waits and launches other jobs.
2. **Calculation Jobs (`calculation`):** These jobs perform the actual calculations.

The `job` section is optional; if it is not defined, the main job will use the same resources specified in the mandatory `calculation` section.

`````{dropdown} Example of cluster section
:color: info
:animate: fade-in-slide-down
:icon: rocket

````{tab} FrontEnd
```yaml
cluster:
  options:
    calculation: None 

```
````
````{tab} SLURM
```yaml
cluster:
  options:
    calculation:
      partition: uds-hub
      time: "2-00:00:00"
      gpus: 1
      mem: '4G'
      constraint: Ryzen_3975WX 
    job:
      partition: uds-hub
      time: "2-00:00:00"
      mem: '1G'
      cpus-per-task: 2

```
````
`````

## `extra_directives` (optional)

### `dependencies` (optional)

This is a list of executable commands that should be run before any `gmx` command. These commands inject the necessary dependencies to ensure the proper execution of the `gmx` command.

````{dropdown} Example of dependencies section
:color: info
:animate: fade-in-slide-down
:icon: rocket
```yaml
extra_directives:
  dependencies:
    - source /groups/CBG/opt/spack-0.18.1/shared.bash
    - module load gromacs/2022.4
    - module load nvidia/latest
    - export GMX_MAXBACKUP=-1

```
````

### `mdrun` (optional)

The user can customize the `gmx mdrun` command. By default, the command `gmx mdrun -nt {threads} -deffnm {simulation_step_name}` is built. You can adjust the `mdrun` options for the ligand, complex, or both simultaneously using the keywords `ligand`, `complex`, and `all`, respectively.

``````{dropdown} Examples of mdrun section
:color: info
:animate: fade-in-slide-down
:icon: rocket

`````{tab} gmx mdrun [...] -cpi -stepout 5000 -v

Example of full command:

```bash
gmx mdrun -nt 12 -deffnm 01_nvt -cpi -stepout 5000 -v
```

````{tab} Configuration acting only on ligand simulations

```yaml
extra_directives:
  mdrun:
    ligand:
      cpi: True
      stepout: 5000
      v: True
```
````

````{tab} Configuration acting only on complex simulations

```yaml
extra_directives:
  mdrun:
    complex:
      cpi: True
      stepout: 5000
      v: True
```
````

````{tab} Configuration acting on both ligand and complex simulations (1)

```yaml
extra_directives:
  mdrun:
    all:
      cpi: True
      stepout: 5000
      v: True
```
````

````{tab} Configuration acting on both ligand and complex simulations (2)

```{tip}
Note that the options for the ligand and complex are independent, meaning that you can set different options for the ligand and complex simulations.
```

```yaml
extra_directives:
  mdrun:
    ligand:
      cpi: True
      stepout: 5000
      v: True
    complex:
      cpi: True
      stepout: 5000
      v: True
```
````

`````


`````{tab} gmx mdrun [...] -ntmpi 1

Example of full command

```bash
gmx mdrun -nt 12 -deffnm 01_nvt -ntmpi 1
```

````{tab} Configuration acting only on ligand simulations

```yaml
extra_directives:
  mdrun:
    ligand:
      ntmpi: 1
```
````

````{tab} Configuration acting only on complex simulations

```yaml
extra_directives:
  mdrun:
    complex:
      ntmpi: 1
```
````

````{tab} Configuration acting on both ligand and complex simulations (1)

```yaml
extra_directives:
  mdrun:
    all:
      ntmpi: 1
```
````

````{tab} Configuration acting on both ligand and complex simulations (2)

```{tip}
Note that the options for the ligand and complex are independent, meaning that you can set different options for the ligand and complex simulations.
```

```yaml
extra_directives:
  mdrun:
    ligand:
      ntmpi: 1
    complex:
      ntmpi: 1
```
````
`````
``````

## `nwindows` (optional)

Number of windows for each step of the perturbation simulations for both the ligand in the solvent and the (membrane) protein-ligand complex. The following is the thermodynamic cycle followed in BindFlow:

```{mermaid}
  graph TB

    no_interacted_lig_in_prot(Non Interacted Restrained Ligand)
    no_coul_lig_in_prot(Restrained Ligand without Coulomb)
    lig_in_prot(Restrained Interacted Ligand)
    free_lig_in_prot(Free Fully Interacted Ligand)
    
    free_lig_in_water(Free Fully Interacted Ligand)
    free_no_coul_lig_in_water(Free Ligand without Coulomb)
    free_no_interacted_lig_in_water(Free Non Interacted Ligand)
    no_interacted_lig_in_water(Non Interacted Restrained Ligand)

  subgraph In the Protein Pocket - 'protein'
    direction BT
    no_interacted_lig_in_prot -- Activate van der Waal - 'vdw' --> no_coul_lig_in_prot -- Activate Coulomb - 'coul' --> lig_in_prot -- Remove Restraints - 'bonded' --> free_lig_in_prot
  end

  subgraph In Water - 'ligand'
    direction TB
    free_lig_in_water -- Remove Coulomb - 'coul' --> free_no_coul_lig_in_water -- Remove van der Waal - 'vdw' --> free_no_interacted_lig_in_water -- Activate Restraints --> no_interacted_lig_in_water
  end
  
  free_lig_in_water -- dG_binding --> free_lig_in_prot
  no_interacted_lig_in_water --  dG = 0 --> no_interacted_lig_in_prot
```

````{dropdown} Example of nwindows section
:color: info
:animate: fade-in-slide-down
:icon: rocket

These are the default options used internally by BindFlow

```yaml   
nwindows:
  ligand:
      vdw: 11
      coul: 11
  complex:
      vdw: 21
      coul: 11
      bonded: 11
```
````

## `mmpbsa` (optional)

This section is used to set the MM(PB/GB)SA calculations. These parameters are passed to [gmx_MMPBSA](https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/) package through the `.in` file.

````{dropdown} Example of mmpbsa section
:color: info
:animate: fade-in-slide-down
:icon: rocket

In this example, the C2, QH, and IE methods are used to estimate entropy, while the PB and GB methods are applied to calculate the polar component of the solvation free energy. For your application, we recommend benchmarking all possible combinations on a small subset of your data to determine the most suitable methods. Once identified, you can streamline your workflow by using a single method (either PB or GB) for solvation free energy calculations and a single method for entropy estimation.

In certain cases, especially within the same family of compounds and a common receptor where relative comparisons are the primary goal, entropy calculations may not be necessary. Omitting them can significantly improve the efficiency of your production workflow.

```yaml   
mmpbsa:
  general: 
    c2_entropy: 1
    qh_entropy: 1
    interaction_entropy: 1
  pb: {} # enable MMPBSA computation
  gb: {} # enable MMGBSA computation
```
````

## `mdp` (optional)

This section is used to control all Molecular Dynamic Parameters for every single simulation

`````{dropdown} Example of mdp section
:color: info
:animate: fade-in-slide-down
:icon: rocket

Here we are only changing the `nsteps` parameter of some of the involved steps.

````{tab} Equilibration
```yaml   
mdp:
  ligand:
    equi:
      01_nvt:
        nsteps: 25
      prod:
        nsteps: 250
  complex:
    equi:
      04_npt:
        nsteps: 25
      prod:
        nsteps: 250

```
````

````{tab} FEP
```yaml   
mdp:
  ligand:
    fep:
      vdw:
        01_nvt:
          nsteps: 25
        prod:
          nsteps: 250
      coul:
        01_nvt:
          nsteps: 25
        prod:
          nsteps: 250
  complex:
    equi:
    fep:
      vdw:
        02_npt:
          nsteps: 25
        03_npt_norest:
          nsteps: 25
        prod:
          nsteps: 250
      coul:
        02_npt:
          nsteps: 25
        03_npt_norest:
          nsteps: 25
        prod:
          nsteps: 250
      bonded:
        02_npt:
          nsteps: 25
        03_npt_norest:
          nsteps: 25
        prod:
          nsteps: 250

```
````

````{tab} MM(PB/GB)SA
```yaml   
mdp:
  complex:
    mmpbsa:
      prod:
        nsteps: 400

```
````

`````

You can explore what are the steps involved in your calculation:

````{tab} Equilibration steps for membrane protein-ligand system

```python
from bindflow.utils.tools import list_if_file
from bindflow.mdp._path_handler import _TemplatePath
print(list_if_file(_TemplatePath.complex.membrane.equi))
```
````

````{tab} FEP steps for soluble protein-ligand system

```python
from bindflow.utils.tools import list_if_file
from bindflow.mdp._path_handler import _TemplatePath
print(list_if_file(_TemplatePath.complex.soluble.fep))
```
````

````{tab} MM(PB/GB)SA step for membrane protein-ligand system

```python
from bindflow.utils.tools import list_if_file
from bindflow.mdp._path_handler import _TemplatePath
print(list_if_file(_TemplatePath.complex.membrane.mmpbsa))
```
````

````{tab} Equilibration steps for the ligand in water

```python
from bindflow.utils.tools import list_if_file
from bindflow.mdp._path_handler import _TemplatePath
print(list_if_file(_TemplatePath.ligand.equi))
```
````

You can also take a look a the default parameters of the step. In the following example, you can print the parameters for the `prod` step of the membrane protein-ligand complex equilibration phase in the MDP format.

```python
from bindflow.mdp._path_handler import _TemplatePath
from bindflow.mdp.mdp import MDP
print(MDP().from_file(_TemplatePath.complex.membrane.equi + "/prod.mdp").to_string())
```

## Suggested options for MM(PB/GB)SA calculations

The default MDP options are optimized for FEP calculations. However, for MM(PB/GB)SA calculations, we recommend using a less resource-intensive scheme. This approach has been shown to be effective, as demonstrated in the main BindFlow publication.

`````{dropdown} Suggested mmpbsa scheme
:color: info
:animate: fade-in-slide-down
:icon: rocket

Note that we are collecting 20 samples (`samples: 20`), and in the `mdp/complex/prod` step, exactly 20 frames are output (`nsteps / nstxout-compressed`) on each XTC file. This configuration ensures the final workflow requires minimal storage, as its size primarily depends on the XTC files.

Additionally, note that a high amount of memory is allocated for the execution of calculations (`cluster/options/mem = 10GB`). This is necessary because [gmx_MMPBSA](https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/) is memory-intensive, as it uses MPI to process frames in parallel. As a general rule, allocate at least 1 GB of memory per thread specified in the `{py:func}`bindflow.runners.calculate` function.

````{tab} Soluble protein-ligand system
```yaml   
cluster:
  options:
    calculation:
      mem: 10G
samples: 20
mdp:
  complex:
    equi:
      00_min:
        nsteps: 100000
      01_nvt:
          dt: 0.002
          nsteps: 5000
      02_nvt:
          dt: 0.003
          nsteps: 5000
      03_npt:
        dt: 0.003
        nsteps: 7500
      04_npt:
        dt: 0.004
        nsteps: 30000
      prod:
        dt: 0.004
        nsteps: 237500
        nstxout-compressed: 11875
    mmpbsa:
      prod:
        dt: 0.004
        nsteps: 25000
        nstxout-compressed: 1250
mmpbsa:
  general: 
    c2_entropy: 1
  gb: {}
```
````
````{tab} Membrane protein-ligand system
```yaml   
cluster:
  options:
    calculation:
      mem: 10G
samples: 20
mdp:
  complex:
    equi:
      00_min:
        nsteps: 100000
      01_nvt:
          dt: 0.001
          nsteps: 5000
      02_nvt:
          dt: 0.001
          nsteps: 5000
      03_npt:
        dt: 0.001
        nsteps: 5000
      04_npt:
        dt: 0.002
        nsteps: 7500
      05_npt:
        dt: 0.002
        nsteps: 7500
      06_npt:
        dt: 0.003
        nsteps: 15000
      prod:
        dt: 0.004
        nsteps: 237500
        nstxout-compressed: 11875
    mmpbsa:
      prod:
        dt: 0.004
        nsteps: 25000
        nstxout-compressed: 1250
mmpbsa:
  general: 
    c2_entropy: 1
  gb: {}
```
````
`````
