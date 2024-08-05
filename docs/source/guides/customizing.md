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
      <mdrun_keywords_for_ligand_simulation>
    all:
      <mdrun_keywords_for_ligand_and_complex_simulation>      
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

## `cluster`

This section specifies the computational resources required by the selected scheduler (see the [Deploy](#bindflow-deploy) section). It allows you to control the allocated resources for two types of jobs:

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

This is a list of executable commands to be executed before any `gmx` command. It helps inject dependencies required for the proper execution of the `gmx` command.

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