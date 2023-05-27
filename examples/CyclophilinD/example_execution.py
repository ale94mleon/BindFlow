#!/usr/bin/env python3

"""
This script is an example execution running for membrane systems.
"""
import yaml
import glob
from abfe import calculate_abfe

ligand_mols = glob.glob("inputs/ligands/*mol")

with open("config.yml", "r") as c:
    global_config = yaml.safe_load(c)

calculate_abfe(
    protein_pdb_path='inputs/protein.pdb',
    ligand_mol_paths=ligand_mols,
    out_root_folder_path="abfe",
    cofactor_mol_path = 'inputs/dummy_cofactor_30.mol',
    membrane_pdb_path = None,
    hmr_factor = 2.5,
    threads = 12,
    ligand_jobs = 1,
    jobs_per_ligand_job = 100000, # I think tht it is better jobs_per_ligand_job
    replicas = 3,
    submit= True,
    global_config = global_config)


# rm -r abfe/slurm_logs/* abfe/*/*/slurm_logs/* abfe/.snakemake/ abfe/*/*/.snakemake abfe/*/*/ligand/ abfe/*/*/complex/