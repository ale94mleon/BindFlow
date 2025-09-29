#!/usr/bin/env python3
import yaml
import glob
from bindflow.runners import calculate
from bindflow.orchestration.generate_scheduler import SlurmScheduler

calculation_type = 'mmpbsa'

force_fields = {
    "espaloma-0.3.1": {
        "type": "espaloma",
        "code": "espaloma-0.3.1"
    },
}

ligand_files = glob.glob("inputs/guests2/*.sdf")

with open(f"config-{calculation_type}.yml", "r") as c:
    global_config = yaml.safe_load(c)

for ff_id, info in force_fields.items():
    ligands = []
    for ligand_file in ligand_files:
        ligands.append({
            'conf': ligand_file,
            'ff': {
                'type': info["type"],
                'code': info["code"]
            }
        })

    calculate(
        calculation_type=calculation_type,
        protein={
            'conf': 'inputs/host/espaloma-0.3.1/OA.gro',
            'top': 'inputs/host/espaloma-0.3.1/OA.top',
        },
        ligands=ligands,
        out_root_folder_path=f"{calculation_type}2/{ff_id}",
        cofactor=None,
        cofactor_on_protein=True,
        host_name='HOST',
        host_selection='resname HOST',
        membrane=None,
        water_model='amber/tip3p',
        hmr_factor=2.5,
        threads=10,
        num_jobs=100000,
        replicas=3,
        submit=False,
        debug=False,
        job_prefix=f'sampl6-oa.{ff_id}',
        scheduler_class=SlurmScheduler,
        global_config=global_config)
