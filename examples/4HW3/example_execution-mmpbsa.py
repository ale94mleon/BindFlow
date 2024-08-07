#!/usr/bin/env python3
import glob

import yaml

from bindflow.orchestration.generate_scheduler import FrontEnd
from bindflow.runners import calculate

ligand_mols = glob.glob("inputs/ligands/*mol")

with open("config-mmpbsa.yml", "r") as c:
    global_config = yaml.safe_load(c)

calculate(
    calculation_type='mmpbsa',
    protein='inputs/protein.pdb',
    ligands=ligand_mols,
    hmr_factor=2.5,
    threads=4,
    num_jobs=12,
    replicas=1,
    scheduler_class=FrontEnd,
    out_root_folder_path="mmpbsa",
    submit=True,
    global_config=global_config)