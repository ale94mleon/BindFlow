#!/usr/bin/env python3
import glob
import yaml
from bindflow.orchestration.generate_scheduler import FrontEnd
from bindflow.run_mmpbsa import calculate_mmpbsa

ligand_mols = glob.glob("inputs/ligands/*mol")

with open("config-mmpbsa.yml", "r") as c:
    global_config = yaml.safe_load(c)

global_config['extra_directives']['mdrun']['all']['ntmpi'] = 1

calculate_mmpbsa(
    protein='inputs/protein.pdb',
    ligands=ligand_mols,
    out_root_folder_path="mmpbsa",
    threads=4,
    num_jobs=12,
    replicas=1,
    submit=True,
    scheduler_class=FrontEnd,
    global_config=global_config)
