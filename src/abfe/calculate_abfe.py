import glob
import os
from typing import List
from abfe.utils import tools
import copy

# TODO: Ligand and RuleThemAll is using CPUs  and they are only waiting
# THink in a way to connect RuleThemAll to Ligand to avoid the use of its CPUs on
# each Ligand simulation
# Think also about in reduce disk space of the simulation
# Maybe at the end of the simulation one command that tar the files
# Do not export so many frames in the xtc file, their are not needed for the analysis.
# For sure not during equilibration phase, keep a realitive small number of frames
from abfe.orchestration.flow_builder import ligand_flows, approach_flow
from abfe.free_energy import gather_results
def calculate_abfe(
        protein_pdb_path: str,
        ligand_mol_paths: List[str],
        out_root_folder_path: str,
        cofactor_mol_path: str = None,
        cofactor_on_protein:bool = True,
        membrane_pdb_path: str = None,
        hmr_factor: float = 3.0,
        dt_max:float = 0.004, # The maximum integration time in ps for all the steps in the workflow. This will be overwrite by the definitions in the global_config
        threads: int = 8, # This is the maximum number of threads to use on the rules, for example to run gmx mdrun
        ligand_jobs: int = None,# By defaults it will take number of ligands * replicas
        jobs_per_ligand_job: int = 10000, # On each ligand, how many jobs should run in parallel
        replicas: int = 3,
        submit: bool = False,
        debug:bool = False,
        job_prefix = None,
        global_config: dict = {}
        ):
    orig_dir = os.getcwd()

    # Make internal copy of configuration
    _global_config =copy.deepcopy(global_config)
     # Check the validity of the provided user configuration file
    check_config = tools.config_validator(global_config=_global_config)
    if not check_config[0]:
        raise ValueError(check_config[1])
    if hmr_factor > 3:
        raise ValueError(f"{hmr_factor =}. Must be lower than 3 to avoid instabilities")
    elif hmr_factor < 2:
        if dt_max > 0.002:
            raise ValueError(f"{hmr_factor =} and {dt_max =}. Is not compatible. for hmr_factor < 2; dt_max must be <= 0.002 ps")
    # IO:
    # Initialize inputs on config
    _global_config["inputs"] = {}
    _global_config["inputs"]["protein_pdb_path"] = os.path.abspath(protein_pdb_path)
    _global_config["inputs"]["ligand_mol_paths"] = [os.path.abspath(ligand_mol_path) for ligand_mol_path in ligand_mol_paths]

    if not _global_config["inputs"]["ligand_mol_paths"]:
        raise ValueError(f'There were not any ligands or they are not accessible on: {ligand_mol_paths}')

    if cofactor_mol_path:
        _global_config["inputs"]["cofactor_mol_path"] = os.path.abspath(cofactor_mol_path)
    else:
        _global_config["inputs"]["cofactor_mol_path"] = None
    _global_config["cofactor_on_protein"] = cofactor_on_protein
    if membrane_pdb_path:
        _global_config["inputs"]["membrane_pdb_path"] = os.path.abspath(membrane_pdb_path)
    else:
        _global_config["inputs"]["membrane_pdb_path"] = None

    _global_config["hmr_factor"] = hmr_factor
    _global_config["dt_max"] = dt_max
    
    out_root_folder_path = os.path.abspath(out_root_folder_path)
    _global_config["out_approach_path"] = out_root_folder_path

    if job_prefix:
        _global_config["job_prefix"] = f"{job_prefix}."
    else:
        _global_config["job_prefix"] = ""
    
    # This will only be needed for developing propose.
    os.environ['abfe_debug'] = str(debug)

    # Generate output folders
    if not os.path.isdir(_global_config["out_approach_path"]):
        tools.makedirs(_global_config["out_approach_path"])

    # Prepare Input / Parametrize
    os.chdir(_global_config["out_approach_path"])

    _global_config["ligand_names"] = [os.path.splitext(os.path.basename(mol))[0] for mol in _global_config["inputs"]["ligand_mol_paths"]]
    _global_config["ligand_jobs"] = ligand_jobs if (ligand_jobs is not None) else len(_global_config["ligand_names"]) * replicas
    _global_config["jobs_per_ligand_job"] = jobs_per_ligand_job
    _global_config["replicas"] = replicas
    _global_config["threads"] = threads

    print("Prepare")
    print("\tstarting preparing ABFE-ligand file structure")

    ligand_flows(_global_config)

    print("\tStarting preparing ABFE-Approach file structure: ", out_root_folder_path)
    expected_out_paths = int(replicas) * len(_global_config["ligand_names"])

    result_paths = glob.glob(_global_config["out_approach_path"] + "/*/*/dG*csv")

    # Only if there is something missing
    if (len(result_paths) != expected_out_paths):
        print("\tBuild approach struct")
        job_id = approach_flow(global_config=_global_config, submit=submit,)
    else:
        job_id = None
    print("Do")
    print("\tSubmit Job - ID: ", job_id)
    # Final gathering
    print("\tAlready got results?: " + str(len(result_paths)))
    if (len(result_paths) > 0):
        print("Trying to gather ready results", out_root_folder_path)
        gather_results.get_all_dgs(root_folder_path=out_root_folder_path, out_csv=os.path.join(out_root_folder_path, 'abfe_partial_results.csv'))
    os.chdir(orig_dir)
