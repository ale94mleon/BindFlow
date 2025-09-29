#!/usr/bin/env python3
import yaml
import glob
from bindflow.runners import calculate
from bindflow.orchestration.generate_scheduler import SlurmScheduler

calculation_type = 'mmpbsa'

force_fields = {
    "openff_unconstrained-2.0.0": {
        "type": "openff",
        "code": "openff_unconstrained-2.0.0.offxml",
    },
    "espaloma-0.3.1": {
        "type": "espaloma",
        "code": "espaloma-0.3.1",
    },
    "gaff-2.11": {
        "type": "gaff",
        "code": "gaff-2.11",
    },
}

ligand_files = glob.glob("inputs/sdf_split/*.mol")

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
            'conf': 'inputs/protein-amber14-all/protein.gro',
            'top': 'inputs/protein-amber14-all/protein.top',
        },
        ligands=ligands,
        membrane='inputs/membrane.pdb',
        cofactor={
            'conf': 'inputs/cofactor.gro',
            'top': 'inputs/cofactor.top',
            'is_water': False,
        },
        cofactor_on_protein=True,
        water_model='amber/tip3p',
        hmr_factor=2.5,
        out_root_folder_path=f"{calculation_type}/{ff_id}",
        threads=10,
        num_jobs=100000,
        replicas=3,
        scheduler_class=SlurmScheduler,
        debug=False,
        job_prefix=f'a2a.{ff_id}',
        submit=False,
        global_config=global_config)
