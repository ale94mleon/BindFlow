import os
import tarfile
from typing import List, Union
import numpy as np
from abfe.orchestration import generate_snake, generate_scheduler
import json

PathLike = Union[os.PathLike, str, bytes]

def update_nwindows_config(config:dict) -> dict:
    """A simple function to update the config file for the entrance nwindows

    Parameters
    ----------
    config : dict
        The configuration file with or without the nwindows keyword. In case it is present, must be in the shape of:
        'nwindows':{
            'ligand':{
                'vdw': <int>[11],
                'coul': <int>[11],
                },
            'complex':{
                'vdw': <int>[21],
                'coul': <int>[11],
                'bonded': <int>[11]
                }, 
            }

    Returns
    -------
    dict
        The updated config
    """
    nwindows_default = {
        'ligand': {
            'vdw': 11,
            'coul': 11,
        },
        'complex':{
            'vdw': 21,
            'coul': 11,
            'bonded':11,
        },
    }
    if 'nwindows' in config:
        nwindows = config['nwindows']
        for key in ['ligand', 'complex']:
            if key in nwindows:
                nwindows_default[key].update(nwindows[key])
    
    config['nwindows'] = nwindows_default
    return config



def ligand_flows(global_config:dict):
    """Prepare all the directories and files

    Parameters
    ----------
    global_config : dict
        This is the global configuration file. It must have:
        num_max_thread: int, The maximum number of threads to be used on each simulation.
        mdrun: dict: A dict of mdrun keywords to add to gmx mdrun, flag must be passed with boolean values. E.g {'cpi': True} 
        extra_dependencies: A list of dependencies that must be run before gmx mdrun. Useful to launch modules as spack or conda.
        out_approach_path: PathLike: The path where the ligands, their replicas and main configuration files will be written.
        num_jobs: int: Maximum number of jobs to run in parallel (TODO, on the ligands??)
        inputs: {
            'ligand_mol_paths': List[PathLike],
            'protein_pdb_path': PathLike,
            'cofactor_mol_path': PathLike,
            'membrane_pdb_path': PathLike
        }

        Usually this dictionary is created by first invoking :meth:`abfe.calculate_abfe.calculate_abfe`
    """

    # Update (or set) nwindows on global_config.
    global_config = update_nwindows_config(global_config)
    # With this implementation the user can select the number of windows setting them up on the global configuration.
    ligand_config = {
        'run_path': None,
        'run_num': None,
        'threads': global_config["threads"],
        'extra_directives': global_config['extra_directives'],
        'num_retries':3,
        'lambdas': {
            'ligand': {
                'vdw': list(np.round(np.linspace(0, 1, global_config['nwindows']['ligand']['vdw']), 2)),
                'coul': list(np.round(np.linspace(0, 1, global_config['nwindows']['ligand']['coul']), 2)),
            },
            'complex': {
                'vdw': list(np.round(np.linspace(0, 1, global_config['nwindows']['complex']['vdw']), 2)),
                'coul': list(np.round(np.linspace(0, 1, global_config['nwindows']['complex']['coul']), 2)),
                'bonded': list(np.round(np.linspace(0, 1, global_config['nwindows']['complex']['bonded']), 2)),
            },
        },
    }
    try:
        ligand_config['mdp'] = global_config['mdp']
    except KeyError:
        pass

    for input_ligand_path in global_config["inputs"]["ligand_mol_paths"]:
        ligand_name = os.path.splitext(os.path.basename(input_ligand_path))[0]
        out_ligand_path = os.path.join(global_config["out_approach_path"],  str(ligand_name))

        # Make directories on demand
        if (not os.path.exists(out_ligand_path)):
            os.mkdir((out_ligand_path))

        out_ligand_input_path = os.path.join(out_ligand_path, "input")
        if (not os.path.isdir(out_ligand_input_path)):
            os.mkdir(out_ligand_input_path)
        
        # Archive original files      
        with tarfile.open(os.path.join(out_ligand_input_path, 'orig_in.tar.gz'), "w:gz") as tar:
            tar.add(input_ligand_path, arcname=os.path.basename(input_ligand_path))
            tar.add(global_config["inputs"]["protein_pdb_path"],arcname=os.path.basename(global_config["inputs"]["protein_pdb_path"]))
            if global_config["inputs"]["cofactor_mol_path"]:
                tar.add(global_config["inputs"]["cofactor_mol_path"],arcname=os.path.basename(global_config["inputs"]["cofactor_mol_path"]))
            if global_config["inputs"]["membrane_pdb_path"]:
                tar.add(global_config["inputs"]["membrane_pdb_path"],arcname=os.path.basename(global_config["inputs"]["membrane_pdb_path"]))
        
        # Update ligand configuration file
        tmp_ligand_config = ligand_config.copy()
        tmp_ligand_config['input_data_path'] = out_ligand_input_path

        # Build the replicas
        for num_replica in range(1, global_config["replicas"] + 1):
            ligand_rep_name = f"{ligand_name}.rep{num_replica}"
            out_replica_path = os.path.join(out_ligand_path, str(num_replica))

            if (not os.path.isdir(out_replica_path)):
                os.mkdir(out_replica_path)

            # set global files:
            snake_path = os .path.join(out_replica_path, "Snakefile")
            conf_path = os .path.join(out_replica_path, "snake_conf.json")

            # Update ligand configuration file
            tmp_ligand_config.update({
                'run_path':out_replica_path,
                'run_num': num_replica,
            })
            # Write the new ligand configuration
            with open(conf_path, "w") as out_IO:
                json.dump(tmp_ligand_config, out_IO, indent=4)
            
            # Generate the snakefile inside each replica
            generate_snake.generate_snake_file(out_file_path=snake_path, conf_file_path=conf_path)

            scheduler = generate_scheduler.create_scheduler(
                scheduler_type = global_config["cluster"]["type"],
                cluster_config = global_config["cluster"]["options"]["calculation"],
                out_dir = out_replica_path,
                prefix_name = ligand_rep_name,
                snake_executor_file = 'job.sh')
            
            scheduler.build_snakemake(jobs = global_config["jobs_per_ligand"])

def approach_flow(global_config: dict, submit=False):
    
    out_path = global_config["out_approach_path"]
    snake_path = out_path + "/Snakefile"
    approach_conf_path = out_path + "/snake_conf.json"

    approach_config = {
        "out_approach_path": global_config["out_approach_path"],
        "inputs": global_config["inputs"],
        "ligand_names": global_config["ligand_names"],
        "replicas": global_config["replicas"],
        'threads': global_config["threads"],
        "hmr_factor": global_config["hmr_factor"]
    }
    with open(approach_conf_path, "w") as out_IO:
        json.dump(approach_config, out_IO, indent=4)

    generate_snake.generate_approach_snake_file(out_file_path=snake_path,
                                                conf_file_path=approach_conf_path)

    # TODO this LIgand si also mostly waiting not performing heavy calculations. Think on a better strategy,
    # I think that in this case we should use the resources of job, but I have to take a look
    scheduler = generate_scheduler.create_scheduler(
        scheduler_type = global_config["cluster"]["type"],
        cluster_config = global_config["cluster"]["options"]["calculation"],
        out_dir = out_path,
        prefix_name = "Ligand",
        snake_executor_file = 'job.sh')
    
    scheduler.build_snakemake(jobs = global_config["ligand_jobs"])

    # Check for extra definitions
    if 'job' in global_config["cluster"]["options"]:
        job_cluster_config = global_config["cluster"]["options"]["job"]
    else:
        job_cluster_config = None
    if submit:
        only_build = False
    else:
        only_build = True
    job_id = scheduler.submit(new_cluster_config = job_cluster_config, only_build = only_build)
    return job_id
