#!/usr/bin/env python3
import glob
import yaml
from bindflow.runners import calculate

ligand_mols = glob.glob("inputs/ligands/*mol")

with open("config.yml", "r") as c:
    global_config = yaml.safe_load(c)

calculate(
    calculation_type='fep',
    protein='inputs/protein.pdb',
    ligands=ligand_mols,
    out_root_folder_path="fep",
    hmr_factor=2.5,
    threads=12,
    num_jobs=100000,
    replicas=3,
    submit=True,
    global_config=global_config)
