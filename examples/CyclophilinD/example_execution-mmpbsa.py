#!/usr/bin/env python3

"""
This script is an example execution running for membrane systems.
"""
import yaml
import glob
from bindflow.run_mmpbsa import calculate_mmpbsa

ligand_mols = glob.glob("inputs/ligands/*mol")

with open("config-mmpbsa.yml", "r") as c:
    global_config = yaml.safe_load(c)

global_config['extra_directives']['mdrun']['all']['ntmpi'] = 1

calculate_mmpbsa(
    protein='inputs/protein.pdb',
    ligands=ligand_mols,
    out_root_folder_path="mmpbsa",
    cofactor='inputs/dummy_cofactor_30.mol',
    membrane=None,
    hmr_factor=2.5,
    threads=12,
    num_jobs=100000,
    replicas=1,
    samples=2,
    submit=False,
    global_config=global_config)


# rm -r abfe/slurm_logs/* abfe/*/*/slurm_logs/* abfe/.snakemake/ abfe/*/*/.snakemake abfe/*/*/ligand/ abfe/*/*/complex/