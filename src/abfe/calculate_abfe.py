import glob
import os
from typing import List
from abfe.utils.tools import config_validator

from abfe.orchestration.flow_builder import ligand_flows, approach_flow
from abfe.scripts import final_receptor_results
def calculate_abfe(
        protein_pdb_path: str,
        ligand_mol_paths: List[str],
        out_root_folder_path: str,
        approach_name: str = "",
        cofactor_mol_path: str = None,
        membrane_pdb_path: str = None,
        hmr_factor: float = 3.0,
        threads: int = 8, # This is the maximum number of threads to use on the rules, for example to run gmx mdrun
        ligand_jobs: int = None,# By defaults it will take number of ligands * replicas
        jobs_per_ligand_job: int = 10000, # On each ligand, how many jobs should run in parallel
        replicas: int = 3,
        submit: bool = False,
        global_config: dict = {}):
    orig_dir = os.getcwd()

    # Check the validity of the provided user configuration file
    check_config = config_validator(global_config=global_config)
    if not check_config[0]:
        raise ValueError(check_config[1])
    if hmr_factor < 2:
        raise ValueError(f'hmr_factor must be equal or higher than 2 (provided {hmr_factor}) to avoid instability during MD simulations. ABFE_workflow uses dt = 4 fs')
    # IO:
    # Initialize inputs on config
    global_config["inputs"] = {}
    global_config["inputs"]["protein_pdb_path"] = os.path.abspath(protein_pdb_path)
    global_config["inputs"]["ligand_mol_paths"] = [os.path.abspath(ligand_mol_path) for ligand_mol_path in ligand_mol_paths]

    if not global_config["inputs"]["ligand_mol_paths"]:
        raise ValueError(f'There were not any ligands or they are not accessible on: {ligand_mol_paths}')

    if cofactor_mol_path:
        global_config["inputs"]["cofactor_mol_path"] = os.path.abspath(cofactor_mol_path)
    else:
        global_config["inputs"]["cofactor_mol_path"] = None
    if membrane_pdb_path:
        global_config["inputs"]["membrane_pdb_path"] = os.path.abspath(membrane_pdb_path)
    else:
        global_config["inputs"]["membrane_pdb_path"] = None

    global_config["hmr_factor"] = hmr_factor
    
    global_config["approach_name"] = approach_name
    global_config["out_approach_path"] = os.path.abspath(out_root_folder_path)

    

    ## Generate output folders
    for dir_path in [global_config["out_approach_path"]]:
        if (not os.path.isdir(dir_path)):
            os.mkdir(dir_path)

    # Prepare Input / Parametrize
    os.chdir(global_config["out_approach_path"])

    global_config["ligand_names"] = [os.path.splitext(os.path.basename(mol))[0] for mol in global_config["inputs"]["ligand_mol_paths"]]
    global_config["ligand_jobs"] = ligand_jobs if (ligand_jobs is not None) else len(global_config["ligand_names"]) * replicas
    global_config["jobs_per_ligand_job"] = jobs_per_ligand_job
    global_config["replicas"] = replicas
    global_config["threads"] = threads

    print("Prepare")
    print("\tstarting preparing ABFE-ligand file structure")

    ligand_flows(global_config)

    print("\tStarting preparing ABFE-Approach file structure: ", out_root_folder_path)
    expected_out_paths = int(replicas) * len(global_config["ligand_names"])
    # TODO, check this part
    result_paths = glob.glob(global_config["out_approach_path"] + "/*/*/dG*csv")

    # Only if there is something missing
    if (len(result_paths) != expected_out_paths):
        print("\tBuild approach struct")
        job_id = approach_flow(global_config=global_config, submit=submit,)
    print("Do")
    print("\tSubmit Job - ID: ", job_id)
    # Final gathering
    print("\tAlready got results?: " + str(len(result_paths)))
    if (len(result_paths) > 0):
        print("Trying to gather ready results", out_root_folder_path)
        final_receptor_results.get_final_results(out_dir=out_root_folder_path, in_root_dir=out_root_folder_path)

    print()
    os.chdir(orig_dir)
