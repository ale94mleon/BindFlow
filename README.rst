BindFlow
========

**WARNING The repo is currently under development - Test systems do not work, code might be instabil :)**

A snakemake based workflow for Absolute Binding Free Energy calculations using GMX. The workflow can be scaled on Slurm Queuing systems. The here provided Cyclophilin D Test systems and experimental values originate from:

* `Alibay, I.; Magarkar, A.; Seeliger, D.; Biggin, P. C. Evaluating the use of absolute binding free energy in the fragment optimisation process. Communications Chemistry 2022, 5, 105. <https://doi.org/10.1038/s42004-022-00721-4>`__
* `Grädler, U.; Schwarz, D.; Blaesse, M.; Leuthner, B.; Johnson, T. L.; Bernard, F.; Jiang, X.; Marx, A.; Gilardone, M.; Lemoine, H.; Roche, D.; Jorand-Lebrun, C. Discovery of novel Cyclophilin D inhibitors starting from three dimensional fragments with millimolar potencies. Bioorganic Medicinal Chemistry Letters 2019, 29, 126717. <https://doi.org/10.1016/j.bmcl.2019.126717>`__

Here a visualization of the triggered process:

|workflow|

TODO
----

Include user friendly options to use the force fields from openmm it is much much easier. The list of force filed for openmm is much higher and easy to access.
In th4e data directory of openmm.

.. code-block:: python

  import parmed as chem
  import parmed.unit as u

  from simtk.openmm import app
  from simtk import openmm as mm

  pdb = app.PDBFile('/home/ale/GIT/protein-ligand-benchmark/plbenchmark/BindFlowData/mcl1/01_protein/crd/protein.pdb')
  forcefield = app.ForceField('amber99sbildn.xml', 'tip3p.xml')
  system = forcefield.createSystem(pdb.topology,)# nonbondedMethod=app.PME, nonbondedCutoff=1*u.nanometer)

  struct = chem.openmm.load_topology(pdb.topology, system, xyz = pdb.positions)
  for path in ['conf.gro', 'topol.top']:
      struct.save(path, overwrite = True)

Install
-------

We will use `mamba <https://mamba.readthedocs.io/en/latest/>`__. First, you must download `environment.yml <https://github.com/bigginlab/ABFE_workflow/blob/main/environment.yml>`__.

If you do not have ``mamba`` installed, then:

.. code-block:: bash

  conda install mamba -n base -c conda-forge

.. warning::

  You could try also with `conda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/windows.html>`__ but it could take a while.

.. code-block:: bash

  mamba env create -f environment.yml

Usage
-----

An example usage is provided with the `examples/example_execution.sh <https://github.com/bigginlab/ABFE_workflow/blob/main/examples/example_execution.sh>`__.

Additional script information is provided via ``abfe-run -h``:

Output:

.. code-block:: bash

  usage: abfe-run [-h] -p PROTEIN_PDB_PATH -l LIGAND_MOL_DIR -o OUTPUT_DIR_PATH [-c COFACTOR_MOL_PATH] [-m MEMBRANE_PDB_PATH] [-hmr HMR_FACTOR] [-nr NUMBER_OF_REPLICATES]
                  [-njr NUMBER_OF_PARALLEL_RECEPTOR_JOBS] [-njl NUMBER_OF_PARALLEL_LIGAND_JOBS] [-ncl NUMBER_OF_CPUS_PER_LIGAND_JOB] [-sc SLRUM_CONFIG] [-submit] [-v]

  optional arguments:
    -h, --help            show this help message and exit
    -p PROTEIN_PDB_PATH, --protein_pdb_path PROTEIN_PDB_PATH
                          Input protein pdb file path
    -l LIGAND_MOL_DIR, --ligand_mol_dir LIGAND_MOL_DIR
                          Input ligand(s) mol file path
    -o OUTPUT_DIR_PATH, --output_dir_path OUTPUT_DIR_PATH
                          Output approach folder
    -c COFACTOR_MOL_PATH, --cofactor_mol_path COFACTOR_MOL_PATH
                          Input cofactor(s) mol file path
    -m MEMBRANE_PDB_PATH, --membrane_pdb_path MEMBRANE_PDB_PATH
                          Input membrane pdb file path. WARNING: The CRYST1 information of this PDB will be used for solvating the system.The protein-membrane system MUST be aligned
                          to the Z-axis (needed for pressure couple scheme).CHARMM-GUI is a good option to get this file.
    -hmr HMR_FACTOR, --hmr_factor HMR_FACTOR
                          The Hydrogen Mass Repartition factor to use. 4 fs of integration time step will be used no matter what hmf_factor is provided. Values greater than 2 are
                          advised, if not the system may be unstable.
    -nr NUMBER_OF_REPLICATES, --number_of_replicates NUMBER_OF_REPLICATES
                          Number of replicates
    -njr NUMBER_OF_PARALLEL_RECEPTOR_JOBS, --number_of_parallel_receptor_jobs NUMBER_OF_PARALLEL_RECEPTOR_JOBS
                          Number of jobs in parallel for receptor workflow
    -njl NUMBER_OF_PARALLEL_LIGAND_JOBS, --number_of_parallel_ligand_jobs NUMBER_OF_PARALLEL_LIGAND_JOBS
                          Number of jobs in parallel for ligand workflow
    -ncl NUMBER_OF_CPUS_PER_LIGAND_JOB, --number_of_cpus_per_ligand_job NUMBER_OF_CPUS_PER_LIGAND_JOB
                          Number of cpus per ligand job
    -sc SLURM_CONFIG, --slrum_config SLURM_CONFIG
                          This is the configuration YAML file of your Slurm cluster. If nothing is provided: partition = cpu time=60:00:00 mem=5000
    -submit               Will automatically submit the ABFE calculations
    -v, --version         show program's version number and exit

Input
-----

The input is suggested to be structured as follows for the command line option:

::

  inputs
  ├── dummy_cofactor_23.mol
  ├── ligands
  │   ├── inhibitor_11.mol
  │   ├── inhibitor_12.mol
  │   ├── inhibitor_17.mol
  │   ├── inhibitor_24.mol
  │   ├── inhibitor_28.mol
  │   ├── inhibitor_2.mol
  │   ├── inhibitor_3.mol
  │   ├── inhibitor_4.mol
  │   ├── inhibitor_6.mol
  │   ├── inhibitor_9.mol
  │   └── ligand.mol
  ├── membrane.pdb
  └── protein.pdb

Running
-------

If the input is set-up correctly and can be parsed, give it a run! (if you want to do the calculation don't forget to `submit`)

Running an ABFE Campaign from Bash:

.. code-block:: bash

  conda activate abfe
  abfe-run -p <path>/receptor.pdb -l <path>/myligands -o <path>/Out -submit

Running an ABFE Campaign from Python

.. code-block:: python

  import glob
  from abfe import calculate_abfe

  ligand_mols = glob.glob("inputs/ligands/*mol")

  out_folder = "abfe"

  calculate_abfe(
      protein_pdb_path='inputs/protein.pdb',
      ligand_mol_paths=ligand_mols,
      out_root_folder_path="abfe",
      membrane_pdb_path = 'inputs/membrane.pdb',
      cofactor_mol_path = 'inputs/dummy_cofactor_23.mol',
      hmr_factor = 3,
      approach_name = "",
      n_cores_per_job= 8,
      num_jobs_receptor_workflow= None,
      num_jobs_per_ligand= 40,
      num_replicas = 3,
      submit= False,
      cluster_config = {})


..  |workflow|  image:: https://github.com/bigginlab/ABFE_workflow/blob/main/.img/full_snakemake_DAG.png?raw=true
    :target: https://github.com/bigginlab/ABFE_workflow/blob/main/.img/full_snakemake_DAG.png
    :alt: logo