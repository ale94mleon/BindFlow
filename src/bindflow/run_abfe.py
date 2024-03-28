import copy
import glob
import os
from typing import List, Union
from warnings import warn

from bindflow._version import __version__
from bindflow.free_energy import gather_results
# TODO: Ligand and RuleThemAll is using CPUs  and they are only waiting
# THink in a way to connect RuleThemAll to Ligand to avoid the use of its CPUs on
# each Ligand simulation
# Think also about in reduce disk space of the simulation
# Maybe at the end of the simulation one command that tar the files
# Do not export so many frames in the xtc file, their are not needed for the analysis.
# For sure not during equilibration phase, keep a relative small number of frames
from bindflow.orchestration.flow_builder import approach_flow
from bindflow.utils import tools

PathLike = Union[os.PathLike, str, bytes]


def input_helper(arg_name: str, user_input: Union[tools.PathLike, dict, None], default_ff: Union[tools.PathLike, str],
                 default_ff_type: Union[str, None] = None, optional: bool = False) -> dict:
    """This helper function is called inside abfe.calculate_ABFE to check for the inputs: protein, ligands, membrane and cofactor

    Parameters
    ----------
    arg_name : str
        The name of the part of the system. It is just used for to print information in case of error
    user_input : Union[tools.PathLike, dict, None]
        The user input provided
    default_ff : Union[tools.PathLike, str]
        A code of the force field. Internally it will be check if [default_ff].ff exist as a directory. This allow a much bigger flexibility
        on the use of different force fields that do not come with the GROMACS distribution by default
    default_ff_type : Union[tools.PathLike, str]
        This is used for the small molecules. It must be openff, gaff or espaloma (case insensitive).
        If it is provided, default_ff will NOT be used and set to None.
        During the building of the system, it will be converted internally as:
            * openff -> openff_unconstrained-2.0.0.offxml
            * gaff -> gaff-2.11
            * espaloma -> espaloma-0.3.1
    optional : bool, optional
        if the arguments under analysis is optional or not, by default False

    Returns
    -------
    dict
        A dictionary with keywords: conf[configuration file], top[GROMACS topology file],
        ff:code[force field code], path[absolute path in case the directory exists]

    Raises
    ------
    ValueError
        if user_input is None but optional is False
    FileNotFoundError
        The configuration file is not found even when some path was provided
    ValueError
        In case conf is not provided when user_input is a dict and optional is False
    FileNotFoundError
        The configuration file is not found when user_input is suppose to be a path
    """
    valid_ff_types = ['openff', 'gaff', 'espaloma']

    if default_ff_type:
        default_ff_type = str(default_ff_type).lower()
        if default_ff_type not in valid_ff_types:
            raise ValueError(f"{default_ff_type =} is not valid. Choose from {valid_ff_types}")

    if not user_input:
        if optional:
            return None
        else:
            raise ValueError(f"{arg_name =} was set with {user_input =} but {optional =}")
    else:
        internal_dict = {
            'conf': None,
            # This must be a single file topology with all the force field information
            # without positional restraint definition for the heavy atoms, thi will be generated internally.
            'top': None,
            'ff': {
                'code': default_ff,
            }
        }
        if default_ff_type:
            internal_dict['ff']['type'] = default_ff_type
            internal_dict['ff']['code'] = None

        if isinstance(user_input, dict):
            tools.recursive_update_dict(internal_dict, user_input)

            # Convert to absolute paths
            if internal_dict['conf']:
                if not os.path.exists(internal_dict['conf']):
                    raise FileNotFoundError(f"{internal_dict['conf'] = } is not accessible.")
                internal_dict['conf'] = os.path.abspath(internal_dict['conf'])
            else:
                if not optional:
                    raise ValueError(f'conf must be provided on the `{arg_name}` entry when a dictionary is used')

            if internal_dict['top']:
                if not os.path.exists(internal_dict['top']):
                    raise FileNotFoundError(f"{internal_dict['top'] = } is not accessible.")
                internal_dict['top'] = os.path.abspath(internal_dict['top'])

        # This is the case that only a path was provided
        else:
            if not os.path.exists(user_input):
                raise FileNotFoundError(f"On {arg_name} entry; {user_input = } is not accessible")
            internal_dict['conf'] = os.path.abspath(user_input)
        return copy.deepcopy(internal_dict)


def calculate_abfe(
        protein: Union[tools.PathLike, dict],  # conf, top, ff
        ligands: Union[tools.PathLike, List[dict]],
        out_root_folder_path: tools.PathLike,
        # You can also specify the keyword is_water in case that the cofactor is a water system,
        # that will change the settles section to triangular constraints. That is needed for compatibility with GROMACS
        cofactor: Union[tools.PathLike, dict, None] = None,
        # this is to the correct group on the thermostat
        cofactor_on_protein: bool = True,
        membrane: Union[tools.PathLike, dict, None] = None,
        # For provided topologies if hmr_factor is set, it will pass any way.
        # So for topology files with already HMR, this should be None.
        # And all the topologies should be provided
        # protein, cofactors, membrane, ligands with the HMR already done
        hmr_factor: Union[float, None] = 3.0,
        # The water force field to use, by default amber/tip3p.
        # if you would likle to use the flexible definition of the CHARMM TIP3P
        # you must define FLEXIBLE and CHARMM_TIP3P in the define statement of the mdp file
        water_model: str = 'amber/tip3p',
        custom_ff_path: Union[None, PathLike] = None,
        # The maximum integration time in ps for all the steps in the workflow.
        # This will be overwrite by the definitions in the global_config
        dt_max: float = 0.004,
        # This is the maximum number of threads to use on the rules, for example to run gmx mdrun
        threads: int = 8,
        # Maximum number of jobs to run in parallel
        num_jobs: int = 10000,
        replicas: int = 3,
        submit: bool = False,
        debug: bool = False,
        job_prefix: Union[None, str] = None,
        global_config: dict = {}
        ):
    print(f"You are using BindFlow: {__version__}.")
    orig_dir = os.getcwd()

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

    # IO:
    # Initialize inputs on config
    _global_config["calculation_type"] = 'fep'  # To leave room for other type of calculations
    _global_config["inputs"] = {}
    _global_config["inputs"]["protein"] = input_helper(arg_name='protein', user_input=protein, default_ff='amber99sb-ildn', optional=False)
    _global_config["inputs"]["ligands"] = [input_helper(arg_name='ligand', user_input=ligand,
                                                        default_ff=None, default_ff_type='openff', optional=False)
                                           for ligand in ligands]
    _global_config["inputs"]["cofactor"] = input_helper(arg_name='cofactor', user_input=cofactor, default_ff=None,
                                                        default_ff_type='openff', optional=True)
    _global_config["inputs"]["membrane"] = input_helper(arg_name='membrane', user_input=membrane, default_ff='Slipids_2020', optional=True)

    _global_config["cofactor_on_protein"] = cofactor_on_protein
    _global_config["hmr_factor"] = hmr_factor
    _global_config["custom_ff_path"] = custom_ff_path
    # TODO, for now I will hard code this section becasue I am modifying the topology with some parameters for the water in preparation.gmx_topology
    _global_config["water_model"] = water_model
    _global_config["dt_max"] = dt_max

    out_root_folder_path = os.path.abspath(out_root_folder_path)
    _global_config["out_approach_path"] = out_root_folder_path

    if job_prefix:
        _global_config["job_prefix"] = f"{job_prefix}"
    else:
        _global_config["job_prefix"] = ""

    # This will only be needed for developing propose.
    os.environ['abfe_debug'] = str(debug)

    # Generate output folders
    if not os.path.isdir(_global_config["out_approach_path"]):
        tools.makedirs(_global_config["out_approach_path"])

    # Prepare Input / Parametrize
    os.chdir(_global_config["out_approach_path"])

    _global_config["ligand_names"] = [os.path.splitext(os.path.basename(mol['conf']))[0] for mol in _global_config["inputs"]["ligands"]]
    _global_config["num_jobs"] = num_jobs
    _global_config["replicas"] = replicas
    _global_config["threads"] = threads

    print("Prepare")
    print("\tstarting preparing ABFE-ligand file structure")

    print("\tStarting preparing ABFE-Approach file structure: ", out_root_folder_path)
    if not _global_config["ligand_names"]:
        raise ValueError("No ligands found")

    expected_out_paths = int(replicas) * len(_global_config["ligand_names"])

    result_paths = glob.glob(_global_config["out_approach_path"] + "/*/*/dG*csv")

    # Only if there is something missing
    if (len(result_paths) != expected_out_paths):
        print("\tBuild approach struct")
        job_id = approach_flow(global_config=_global_config, submit=submit)
    else:
        job_id = None
    print("Do")
    print("\tSubmit Job - ID: ", job_id)
    # Final gathering
    print("\tAlready got results?: " + str(len(result_paths)))
    if (len(result_paths) > 0):
        print("Trying to gather ready results", out_root_folder_path)
        gather_results.get_all_dgs(root_folder_path=out_root_folder_path, out_csv=os.path.join(out_root_folder_path, 'abfe_partial_results.csv'))
        gather_results.get_raw_data(root_folder_path=out_root_folder_path, out_csv=os.path.join(out_root_folder_path, 'abfe_partial_results_raw.csv'))
    os.chdir(orig_dir)
