import copy
import glob
import os
from pathlib import Path
from typing import List, Union
from warnings import warn

from bindflow._gmx_check import check_gromacs_installation
from bindflow._version import __version__
from bindflow.free_energy import gather_results
from bindflow.orchestration.flow_builder import approach_flow
from bindflow.orchestration.generate_scheduler import Scheduler, SlurmScheduler
from bindflow.utils import tools

import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


PathLike = Union[os.PathLike, str, bytes]


def calculate(
        calculation_type: str,
        protein: Union[tools.PathLike, dict],
        ligands: Union[tools.PathLike, List[dict]],
        membrane: Union[tools.PathLike, dict, None] = None,
        cofactor: Union[tools.PathLike, dict, None] = None,
        cofactor_on_protein: bool = True,
        water_model: str = 'amber/tip3p',
        custom_ff_path: Union[None, PathLike] = None,
        host_name: str = 'Protein',
        host_selection: str = 'protein and name CA',
        fix_protein: bool = True,
        hmr_factor: Union[float, None] = 2.5,
        dt_max: float = 0.004,
        threads: int = 12,
        num_jobs: int = 10000,
        replicas: int = 3,
        scheduler_class: Scheduler = SlurmScheduler,
        debug: bool = False,
        job_prefix: Union[None, str] = None,
        out_root_folder_path: tools.PathLike = 'bindflow-out',
        submit: bool = False,
        global_config: dict = {}
        ) -> None:
    """Main function of BindFlow to execute the workflow

    Parameters
    ----------
    calculation_type : str

        Any of (case-insensitive):

            * "fep": For Free Energy Perturbation simulations
            * "mmpbsa": For Molecular  Molecular Mechanic Poisson-Boltzmann/Generalized-Born Surface Area MM(PB/GB)SA simulations

    protein : Union[tools.PathLike, dict]
        This could be the path to the PDB file of the protein which will be processed through
        GMX with amber99sb-ildn; or a dictionary with the specific definition of the protein.

        In case a dictionary is provided, it should have:

            * conf -> The path of the protein PDB/GRO file [mandatory]

            * top -> GROMACS topology [optional], by default None.
            Should be a single file topology with all the force field
            information and without the position restraint included. However, in case,
            you need to use an include statement such as:

                include "./charmm36-jul2022.ff/forcefield.itp"

            You must change the statement to the absolute path:

                include "{prefix of the absolute path}/charmm36-jul2022.ff/forcefield.itp"

            And copy the charmm36-jul2022.ff to custom_ff_path and set this parameter accordingly. If not
            you may get some errors about files not founded. The force field directory
            must end with the suffix ".ff".

            * ff
                * code -> GMX force field code [optional], by default amber99sb-ildn
                You can use your custom force field, but custom_ff_path must be provided

    ligands : Union[tools.PathLike, List[dict]]
        This is a list of either path to the MOL/SDF file of the ligands which will be processed through
        TOFF with openff_unconstrained-2.0.0.offxml; or a dictionary which expose more options to use with
        the TOFF Python library; or a combination of both.

        In case the element is a dictionary, it should have:

            * conf -> The path of the small molecule MOL/SDF file [mandatory]. In case that top is provided,
            this must be a .gro, a ValueError will be raised if it is not the case
            the molecule will not get its parameters.

            * top -> GROMACS topology [optional]. Must be a single file topology with all the force field
            information and without the position restraint included, by default None

            * ff:

                * type -> openff, gaff or espaloma

                * code -> force field code [optional], by default depending on type

                    * openff -> openff_unconstrained-2.0.0.offxml

                    * gaff -> gaff-2.11

                    * espaloma -> espaloma-0.3.1

                With this parameter you can access different small molecule force fields

    membrane : Union[tools.PathLike, dict, None], optional
        This is either None (default); a path to the PDB file of the membrane which will be processed
        through GMX with SLipid2020; or a dictionary with the specific definition of the protein.

        In case a dictionary is provided, it should have:

            * conf -> The path of the membrane PDB file [mandatory]. If provided, the PDB must have a
            correct definition of the CRYST1. This information will be used for the solvation step.
            The membrane must be already correctly placed around the protein. Servers like CHARM-GUI
            can be used on this step.

            * top -> GROMACS topology [optional], by default None.
            Should be a single file topology with all the force field
            information and without the position restraint included. However, in case,
            you need to use an include statement such as:

                include "./amber-lipids14.ff/forcefield.itp"

            You must change the statement to the absolute path:

                include "{prefix of the absolute path}/amber-lipids14.ff/forcefield.itp"

            And copy theamber-lipids14.ff to custom_ff_path and set this parameter accordingly. If not
            You may get some errors about files not founded. The force field directory
            must end with the suffix ".ff".

            * ff

                * code -> GMX force field code [optional], by default Slipids_2020
                You can use yoru custom force field, but custom_ff_path must be provided

    cofactor : Union[tools.PathLike, dict, None], optional
        This is either None (default); a path to the MOL/SDF file of the ligands which will be processed
        through TOFF with openff_unconstrained-2.0.0.offxml; or a dictionary which expose more options
        to use with the TOFF Python library

        In case the element is a dictionary, it should have:

            * conf -> The path of the small molecule MOL/SDF file [mandatory]. In case that top is provided,
            this must be a .gro, a ValueError will be raised if it is not the case
            the molecule will not get its parameters.

            * top -> GROMACS topology [optional]. Must be a single file topology with all the force field
            information and without the position restraint included, by default None

            * ff:

                * type -> openff, gaff or espaloma

                * code -> force field code [optional], by default depending on type

                    * openff -> openff_unconstrained-2.0.0.offxml

                    * gaff -> gaff-2.11

                    * espaloma -> espaloma-0.3.1

                With this parameter you can access different small molecule force fields

            * is_water -> If presents and set to True; it is assumed that this is a water system
            and that will change the settles section (if any) to tip3p-like triangular constraints.
            This is needed for compatibility with GROMACS. Check here:
            https://gromacs.bioexcel.eu/t/how-to-treat-specific-water-molecules-as-ligand/3470/9

    cofactor_on_protein : bool, optional
        It is used during the index generation for membrane systems. It only works if cofactor_mol is provided.
        If True, the cofactor will be part of the protein and the ligand
        if False will be part of the solvent and ions. This is used mainly for the thermostat. By default True

    water_model : str, optional
        The water force field to use, by default amber/tip3p.
        if you would like to use the flexible definition of the CHARMM TIP3P
        you must define FLEXIBLE and CHARMM_TIP3P in the define statement of the mdp file

    custom_ff_path : Union[None, PathLike], optional
        All the custom force field must be in this directory. The class will set:

            os.environ["GMXLIB"] = os.path.abspath(custom_ff_path)

    host_name : str, optional
        The group name for the host in the configuration file, by default "Protein".
         This is used for making index, solvate the system and working with trajectories

    host_selection : str, optional
        MDAnalysis selection to define the host (receptor or protein), by default 'protein and name CA'.
        This is used for Boresch restraint detection.

    fix_protein : bool, optional
        If True, `pdbfixer` will be applied with flags `--add-atoms=all --replace-nonstandard` and `gmx editconf`
        will the `-ignh` flag. This is needed to avoid possible issues when processing the structure through
        GROMACS. To kept an specific protonation state is advised to input the full definition of
        the protein (.top, .gro) or a PDB with the atom-naming (mainly H-naming) consistent with your selected
        force field. This should be used for protein mainly, by default True

    hmr_factor : Union[float, None], optional
        The Hydrogen Mass Factor to use, by default 2.5.

        .. warning::
            For provided topologies if hmr_factor is set, it will pass any way.
            So for topology files with already HMR, this should be None.
            And all the topologies should be provided
            protein, cofactors, membrane, ligands with the HMR already done

    dt_max : float, optional
        This is the maximum integration time step that will be used by any MD simulation step
        This will be override by the specific MDP step definition through the  the definitions
        in the global_config, by default 0.004
    threads : int, optional
        This is the maximum number of CPUs/threads to use by any Snakemake rule. E.g. `gmx mdrun` will run with this amount of threads, by default 12
    num_jobs : int, optional
        This is the maximum Snakemake concurrent jobs, by default 10000.
        When you launch in a HPC (e.g.Slurm) you can use (if your system allows it) a high number; In this case Snakemake counts as running
        jobs both those ones actually running and the pending ones.
        In the other hand, if (for testing or any other use) the FrontEnd is been used, this parameter should be set to the amount of CPUs that
        you would like to allocate for the entire workflow. This will prevent to overheat your machine.
        For example in a workstation of 12 CPus, if you set threads = 4, then num_jobs should be 3.
    replicas : int, optional
        The number of independent repeats of the entire workflow (the building of the system is not repeated), by default 3

    scheduler_class : Scheduler, optional
        This is a class to schedule the jobs and specify how to handle computational resources, by default SlurmScheduler

        The module :mod:`bindflow.orchestration.generate_scheduler` presents the template class
        :class:`bindflow.orchestration.generate_scheduler.Scheduler` which can be used to create customized Scheduler based on user needs.
        :mod:`bindflow.orchestration.generate_scheduler` also contains the following functional and already tested schedular:

        #. :class:`bindflow.orchestration.generate_scheduler.SlurmScheduler`: To interact with `Slurm <https://slurm.schedmd.com/documentation.html>`_
        #. :class:`bindflow.orchestration.generate_scheduler.FrontEnd`: To execute the workflow in a frontend-like computer. E.g. LAPTOP, workstation, etc.

    debug : bool, optional
        If True more stuff will be printed, by default False
    job_prefix : Union[None, str], optional
        A prefix to identify the jobs in the HPc cluster queue, by default None
    out_root_folder_path : tools.PathLike
        Where the workflow is going to run, by default bindflow-out
    submit : bool, optional
        If True the workflow will woke alive, by default False

    global_config : dict, optional
        The rest of the configuration and fine tunning of the workflow goes here, by default {}

    Raises
    ------
    ValueError
        In case of invalid global_config
    ValueError
        In case the ligand paths are not found
    ValueError
        In case wrong calculation_type
    RuntimeError
        For incompatible GROMACS version
    """

    logging.info(f"✨ You are using BindFlow: {__version__}")

    if calculation_type.lower() not in ['fep', 'mmpbsa']:
        raise ValueError(f"calculation_type must be one of: [fep, mmpbsa] (case-insensitive).\nProvided: {calculation_type}")
    else:
        calculation_type = calculation_type.lower()

    check_gromacs_installation()
    out_root_folder_path = Path(out_root_folder_path)

    # Make internal copy of configuration
    _global_config = copy.deepcopy(global_config)
    # Check the validity of the provided user configuration file
    check_config = tools.config_validator(global_config=_global_config)
    if not check_config[0]:
        raise ValueError(check_config[1])
    if hmr_factor:
        if hmr_factor > 4:
            warn(f"{hmr_factor =}. It should be lower or equal than 4 (preferred 3) to avoid instabilities")
        elif hmr_factor < 2:
            if dt_max > 0.002:
                warn(f"{hmr_factor =} and {dt_max =}. For hmr_factor < 2; dt_max should be <= 0.002 ps")
    else:
        if dt_max > 0.002:
            warn(f"{hmr_factor =} and {dt_max =}. hmr_factor is not been, therefore dt_max should be <= 0.002 ps")

    # Initialize inputs on config
    _global_config["calculation_type"] = calculation_type
    _global_config["scheduler_class"] = scheduler_class
    _global_config["inputs"] = {}
    _global_config["inputs"]["protein"] = tools.input_helper(arg_name='protein', user_input=protein, default_ff='amber99sb-ildn', optional=False)
    # TODO check that is a list, tuple or string, iterable is nto enough because the dict is an iterable. Not clear how to check for this
    _global_config["inputs"]["ligands"] = [tools.input_helper(arg_name='ligand', user_input=ligand,
                                                              default_ff=None, default_ff_type='openff', optional=False)
                                           for ligand in ligands]
    _global_config["inputs"]["cofactor"] = tools.input_helper(arg_name='cofactor', user_input=cofactor, default_ff=None,
                                                              default_ff_type='openff', optional=True)
    _global_config["inputs"]["membrane"] = tools.input_helper(arg_name='membrane', user_input=membrane, default_ff='Slipids_2020', optional=True)

    _global_config["host_name"] = host_name
    _global_config["host_selection"] = host_selection
    _global_config["fix_protein"] = fix_protein
    _global_config["cofactor_on_protein"] = cofactor_on_protein
    _global_config["hmr_factor"] = hmr_factor
    _global_config["custom_ff_path"] = custom_ff_path
    # TODO, for now I will hard code this section becasue I am modifying the topology with some parameters for the water in preparation.gmx_topology
    _global_config["water_model"] = water_model
    _global_config["dt_max"] = dt_max
    _global_config["out_approach_path"] = os.path.abspath(out_root_folder_path)

    if job_prefix:
        _global_config["job_prefix"] = f"{job_prefix}"
    else:
        _global_config["job_prefix"] = ""

    # This will only be needed for developing propose.
    os.environ['BINDFLOW_DEBUG'] = str(debug)

    # Generate output folders
    if not Path(_global_config["out_approach_path"]).is_dir():
        Path(_global_config["out_approach_path"]).mkdir(exist_ok=True, parents=True)

    # Prepare Input / Parametrize

    _global_config["ligand_names"] = [Path(mol['conf']).stem for mol in _global_config["inputs"]["ligands"]]
    _global_config["num_jobs"] = num_jobs
    _global_config["replicas"] = replicas
    _global_config["threads"] = threads

    # Check default samples for mmpbsa simulations
    if calculation_type == 'mmpbsa':
        if 'samples' in _global_config:
            samples = _global_config['samples']
        else:
            _global_config['samples'] = 20
            samples = 20

    logging.info(f"🏗️  Building file structure for {calculation_type}: {out_root_folder_path}")

    if not _global_config["ligand_names"]:
        raise ValueError("No ligands found")

    if calculation_type == 'fep':
        expected_out_paths = int(replicas) * len(_global_config["ligand_names"])
        result_paths = glob.glob(_global_config["out_approach_path"] + "/*/*/dG*csv")
    elif calculation_type == 'mmpbsa':
        expected_out_paths = replicas * samples * len(_global_config["ligand_names"])
        result_paths = glob.glob(_global_config["out_approach_path"] + "/*/*/complex/mmpbsa/simulation/*/mmxbsa.csv")

    # Only if there is something missing
    if (len(result_paths) != expected_out_paths):
        job_id = approach_flow(global_config=_global_config, submit=submit)
        if job_id:
            logging.info(f"🚀 Submit Job - ID: {job_id}")
        else:
            logging.info("🛰️  BindFlow tasks are not yet submitted")
    else:
        logging.info("✅ All gathering CSV files were generated, nothing to do.")
    if (len(result_paths) > 0):
        print(f"🗃️ Trying to gather {len(result_paths)} ready results on: {out_root_folder_path}")
        if calculation_type == 'fep':
            gather_results.get_all_fep_dgs(root_folder_path=out_root_folder_path,
                                           out_csv=out_root_folder_path/'fep_partial_results.csv')
            gather_results.get_raw_fep_data(root_folder_path=out_root_folder_path,
                                            out_csv=out_root_folder_path/'fep_partial_results_raw.csv')
        elif calculation_type == 'mmpbsa':
            full_df = gather_results.get_raw_mmxbsa_dgs(root_folder_path=out_root_folder_path,
                                                        out_csv=out_root_folder_path/'mmxbsa_partial_results_raw.csv')
            gather_results.get_all_mmxbsa_dgs(full_df=full_df, columns_to_process=None,
                                              out_csv=out_root_folder_path/'mmxbsa_partial_results.csv')
