#!/usr/bin/env python3

"""
This script is an example execution running for membrane systems.
"""
import yaml
import glob
from bindflow import calculate_mmpbsa

ligand_mols = glob.glob("inputs/ligands/*mol")

with open("config-mmpbsa.yml", "r") as c:
    global_config = yaml.safe_load(c)

calculate_mmpbsa(
    protein='inputs/protein.pdb',
    ligands=ligand_mols,
    out_root_folder_path="mmpbsa",
    cofactor=None,
    membrane=None,
    hmr_factor=2.5,
    threads=10,
    num_jobs=10,
    replicas=1,
    samples=20,
    submit=False,
    global_config=global_config)


# rm -r abfe/slurm_logs/* abfe/*/*/slurm_logs/* abfe/.snakemake/ abfe/*/*/.snakemake abfe/*/*/ligand/ abfe/*/*/complex/