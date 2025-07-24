import json
import os
import tarfile
from pathlib import Path
from typing import Union

import numpy as np

from bindflow import rules

PathLike = Union[os.PathLike, str, bytes]


def update_nwindows_config(config: dict) -> dict:
    """A simple function to update the config file for the entrance nwindows

    Parameters
    ----------
    config : dict
        The configuration file with or without the nwindows keyword.
        In case it is present, must be in the shape of:
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
        'complex': {
            'vdw': 21,
            'coul': 11,
            'bonded': 11,
        },
    }
    if 'nwindows' in config:
        nwindows = config['nwindows']
        for key in ['ligand', 'complex']:
            if key in nwindows:
                nwindows_default[key].update(nwindows[key])

    config['nwindows'] = nwindows_default
    return config


def generate_approach_snake_file(out_file_path: str, conf_file_path: str, calculation_type: str) -> None:
    """Used to generate the main Snakefile

    Parameters
    ----------
    out_file_path : str
        Path to write the Snakefile
    conf_file_path : str
        Path of the yml workflow configuration file.
    calculation_type : str
        Either mmpbsa or fep.
    """
    # Sanity check
    valid_calculation_type = ['mmpbsa', 'fep']
    if calculation_type not in valid_calculation_type:
        raise ValueError(f"{calculation_type} is an invalid calculation_type, choose from {valid_calculation_type}")
    file_str = "# Load Config:\n"\
        f"configfile: \'{conf_file_path}\'\n"\
        "from pathlib import Path\n\n"\
        "# Start Flow\n"\
        f"include: \'{rules.super_flow}/Snakefile\'\n\n"\
        "# Specify targets and dependencies\n"\
        "rule RuleThemAll:\n"

    if calculation_type == 'fep':
        file_str += "    input: Path(config[\"out_approach_path\"])/\"fep_results.csv\""
    elif calculation_type == 'mmpbsa':
        file_str += "    input: Path(config[\"out_approach_path\"])/\"mmxbsa_results.csv\""

    with open(out_file_path, 'w') as out:
        out.write(file_str)


def approach_flow(global_config: dict, submit: bool = False) -> str:
    """It controls the rest of the workflows
    that make the actual calculations. It will only hang and wait till the rest
    subprocess finish. In case that cluster/options/job is defined in global_config,
    those options will be used to create the proper cluster submit script, if not
    cluster/option/calculation will be used instead

    Parameters
    ----------
    global_config : dict
        The global configuration. It should contain:
        out_approach_path[PathLike], inputs[dict[dict]], water_model[str],
        host_name[str], host_selection[str] (no needed for mmpbsa),
        cofactor_on_protein[bool], extra_directives[dict], dt_max[float]
        ligand_names[list[str]], replicas[float], threads[int], samples[int] (no needed for fep)
        hmr_factor[float, None], custom_ff_path[str, None], cluster/type[str], cluster/options/calculation[dict]
        num_max_thread: int, The maximum number of threads to be used on each simulation.
        mdrun: dict: A dict of mdrun keywords to add to gmx mdrun, flag must be passed with boolean values. E.g {'cpi': True}
        extra_dependencies: A list of dependencies that must be run before gmx mdrun. Useful to launch modules as spack or conda.
        num_jobs: int: Maximum number of jobs to run in parallel
        cluster/options/job[dict]. The last is optional and will override cluster/options/calculation[dict]
        during submit
    submit : bool, optional
        Submit to the workload manager, by default False

    Returns
    -------
    str
        Some identification of the submitted job. It will depend on how
        the submit method of the corresponded Schedular (:class:`bindflow.orchestration.generate_scheduler.Scheduler`) was implemented
    """
    out_path = Path(global_config["out_approach_path"])
    snake_path = out_path/"Snakefile"
    approach_conf_path = out_path/"snake_conf.json"       

    approach_config = {
        "calculation_type": global_config["calculation_type"],
        "out_approach_path": str(global_config["out_approach_path"]),
        "inputs": global_config["inputs"],
        "water_model": global_config["water_model"],
        "host_name": global_config["host_name"],
        "fix_protein": global_config["fix_protein"],
        "cofactor_on_protein": global_config["cofactor_on_protein"],
        "ligand_names": global_config["ligand_names"],
        "replicas": global_config["replicas"],
        "hmr_factor": global_config["hmr_factor"],
        "custom_ff_path": global_config["custom_ff_path"],
        'threads': global_config['threads'],
        'extra_directives': global_config['extra_directives'],
        'retries': 3,
        'dt_max': global_config['dt_max'],
        # With this implementation the user can select the number of windows setting them up on the global configuration.
    }
    if global_config["calculation_type"] == 'fep':
        # Update number of windows if needed and create the lambda-schedule
        global_config = update_nwindows_config(global_config)
        approach_config['lambdas'] = {
            'ligand': {
                'vdw': list(np.round(np.linspace(0, 1, global_config['nwindows']['ligand']['vdw']), 2)),
                'coul': list(np.round(np.linspace(0, 1, global_config['nwindows']['ligand']['coul']), 2)),
            },
            'complex': {
                'vdw': list(np.round(np.linspace(0, 1, global_config['nwindows']['complex']['vdw']), 2)),
                'coul': list(np.round(np.linspace(0, 1, global_config['nwindows']['complex']['coul']), 2)),
                'bonded': list(np.round(np.linspace(0, 1, global_config['nwindows']['complex']['bonded']), 2)),
            },
        }
        approach_config["host_selection"] = global_config["host_selection"]
    elif global_config["calculation_type"] == 'mmpbsa':
        approach_config["samples"] = global_config["samples"]
        if "mmpbsa" in global_config.keys():
            approach_config["mmpbsa"] = global_config["mmpbsa"]

    # Specify the complex type
    if global_config["inputs"]["membrane"]:
        approach_config["complex_type"] = 'membrane'
    else:
        approach_config["complex_type"] = 'soluble'

    # Add extra mdp options if provided
    try:
        approach_config['mdp'] = global_config['mdp']
    except KeyError:
        pass

    # Just to save the prefix
    if global_config["job_prefix"]:
        approach_config["job_prefix"] = global_config["job_prefix"]

    for ligand_definition in global_config["inputs"]["ligands"]:
        input_ligand_path = Path(ligand_definition['conf'])
        ligand_name = input_ligand_path.stem
        out_ligand_path = Path(global_config["out_approach_path"])/ligand_name

        # Make directories on demand
        out_ligand_path.mkdir(exist_ok=True, parents=True)
        out_ligand_input_path = out_ligand_path/"input"
        out_ligand_input_path.mkdir(exist_ok=True, parents=True)
        (out_ligand_input_path/"complex").mkdir(exist_ok=True, parents=True)
        (out_ligand_input_path/"ligand").mkdir(exist_ok=True, parents=True)

        # Archive original files
        with tarfile.open(out_ligand_input_path/'orig_in.tar.gz', "w:gz") as tar:
            tar.add(input_ligand_path, arcname=input_ligand_path.name)
            tar.add(global_config["inputs"]["protein"]["conf"], arcname=Path(global_config["inputs"]["protein"]["conf"]).name)
            if global_config["inputs"]["cofactor"]:
                tar.add(global_config["inputs"]["cofactor"]["conf"], arcname=Path(global_config["inputs"]["cofactor"]["conf"]).name)
            if global_config["inputs"]["membrane"]:
                tar.add(global_config["inputs"]["membrane"]["conf"], arcname=Path(global_config["inputs"]["membrane"]["conf"]).name)

        # Build the replicas
        for num_replica in range(1, global_config["replicas"] + 1):
            out_replica_path = out_ligand_path/str(num_replica)
            out_replica_path.mkdir(exist_ok=True, parents=True)

    with open(approach_conf_path, "w") as out_IO:
        json.dump(approach_config, out_IO, indent=4)

    generate_approach_snake_file(out_file_path=snake_path, conf_file_path=approach_conf_path, calculation_type=global_config["calculation_type"])

    scheduler_class = global_config['scheduler_class']
    scheduler = scheduler_class(
        # by default, run with the main cluster options
        # only if global_config["cluster"]["options"]["job"] is defined it will change during submit
        cluster_config=global_config["cluster"]["options"]["calculation"],
        out_dir=out_path,
        prefix_name=f"{global_config['job_prefix']}",
        snake_executor_file='job.sh')

    scheduler.build_snakemake(jobs=global_config["num_jobs"])

    # Check for extra definitions
    if 'job' in global_config["cluster"]["options"]:
        job_cluster_config = global_config["cluster"]["options"]["job"]
    else:
        job_cluster_config = None

    # if global_config["cluster"]["options"]["job"] changes during submit the cluster options
    # Execute the pipeline in out_approach_path
    cwd = os.getcwd()
    os.chdir(global_config["out_approach_path"])
    job_id = scheduler.submit(new_cluster_config=job_cluster_config, only_build=not submit, job_prefix=global_config["job_prefix"])
    os.chdir(cwd)
    return job_id
