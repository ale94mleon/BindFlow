import json
from pathlib import Path
from bindflow.utils.tools import List, PathLike, list_if_file

_MDP_PARAM_DEFAULT = {
    "integrator": "steep",
    "emtol": "1000.0",
    "nsteps": "5000",
    "nstlist": "10",
    "cutoff-scheme": "Verlet",
    "rlist": "1.0",
    "vdwtype": "Cut-off",
    "vdw-modifier": "Potential-shift-Verlet",
    "rvdw-switch": "0",
    "rvdw": "1.0",
    "coulombtype": "pme",
    "rcoulomb": "1.0",
    "epsilon-r": "1",
    "epsilon-rf": "1",
    "constraints": "h-bonds",
    "constraint-algorithm": "LINCS"
}


class MDP:
    """Base class to work with MDP files
    """
    def __init__(self, **kwargs):
        self.parameters = dict()
        self._set_default_parameters()
        self.set_parameters(**kwargs)

    def _set_default_parameters(self):
        # Add any default parameters here
        self.parameters = _MDP_PARAM_DEFAULT

    def set_parameters(self, **kwargs):
        kwargs = {key.replace('_', '-'): value for key, value in kwargs.items()}
        self.parameters.update(kwargs)

    def from_file(self, template_filename, clean_current_parameters=True):
        with open(template_filename, 'r') as f:
            lines = f.readlines()
            if clean_current_parameters:
                # Clean all defined parameters
                self.parameters = {}
            for line in lines:
                if line.startswith(';') or line.startswith('#'):
                    continue
                tokens = line.strip().split('=', 1)
                if len(tokens) != 2:
                    continue
                parameter_name = tokens[0].strip().replace('_', '-')
                parameter_value = tokens[1].strip()
                self.parameters[parameter_name] = parameter_value
        return self

    def to_string(self):
        s = ''
        for parameter_name, parameter_value in self.parameters.items():
            s += f'{parameter_name:<40} = {parameter_value}\n'
        return s

    def write(self, filename: str):
        with open(filename, 'w') as f:
            f.write(self.to_string())

    def __repr__(self):
        return f"{self.__class__.__name__}({json.dumps(self.parameters, indent=4)})"


class StepMDP(MDP):
    """This subclass will inherit from :class:`bindflow.mdp.mdp.MDP`
    It is meant to be used in combination with the templates that can
    be access from ``bindflow.mdp.templates.TemplatePath``.
    This class define the method ``set_new_step``. One time initialized,
    the instance could be used to access other steps on the step_path

    Parameters
    ----------
    MDP : :class:`bindflow.mdp.mdp.MDP`
        base MDP class
    """
    def __init__(self, step: str = None, step_path: PathLike = None, **kwargs):
        """Constructor. It is assume a tree directory as:

        .. code-block:: text

            .
            ├── emin.mdp
            ├── npt.mdp
            ├── npt-norest.mpd
            ├── nvt.mdp
            └── prod.mdp

        Parameters
        ----------
        step : str, optional
            the step, basically the name of the mdp file on templates, by default None
        step_path : PathLike, optional
            where to look for the mdp, by default None

        Example
        -------
        .. ipython:: python

            from bindflow.mdp import mdp
            from bindflow.mdp.templates import TemplatePath
            from pathlib import Path
            my_mdp = mdp.StepMDP(step='00_min', step_path=Path(TemplatePath.ligand.fep)/'coul')
            my_mdp.set_parameters(**{"init-lambda-state": "0", "coul-lambdas": "0 0.5 1"})
            print(my_mdp)
            my_mdp.set_new_step(step='01_nvt')
            print(my_mdp.to_string())
        """
        super().__init__(**kwargs)
        self.step = step
        self.step_path = Path(step_path)
        if self.step:
            self.__from_archive()

    def set_new_step(self, step):
        self.__from_archive(explicit_step=step)
        return self

    def __from_archive(self, explicit_step: str = None):
        if explicit_step:
            self.step = explicit_step
        valid_steps = [step.stem for step in list_if_file(self.step_path, ext='.mdp')]
        if self.step not in valid_steps:
            raise ValueError(f"name = {self.step} is not a valid step mdp, must be one of: {valid_steps}")
        self.from_file(self.step_path/f"{self.step}.mdp")


def make_fep_dir_structure(sim_dir: PathLike, template_dir: PathLike, lambda_values: List[float], lambda_type: str, sys_type: str,
                           dt_max: float, mdp_extra_kwargs: dict = None):
    """This function is meant to be used on ``ligand_fep_setup`` and ``complex_fet_setup`` rules.
    It will create the structure of the simulation directory: ``{sim_dir}/simulation/{lambda_type}.{i}/{step}/{step}.mdp``

    Where:

        * i:  init-lambda-state,
        * step: the name of the simulation to carry on

    Parameters
    ----------
    sim_dir : PathLike
        Where the simulation suppose to run
    template_dir : PathLike
        This is the directory that storage the mdp templates: bindflow.mdp.templates.TemplatePath.ligand.fep or
        bindlfow.mdp.templates.TemplatePath.complex.fep
    lambda_values : List[float]
        This is a the list of lambda values to be used inside the mdp on the entrance {lambda_type}-lambdas
    lambda_type : str
        Must be one of the following strings "vdw", "coul", "bonded" (the last is for restraints)
    sys_type : str
        Must one of the following strings "ligand" or "complex". This is used in order to turn on the bonded
        lambdas for the complex simulations
    mdp_extra_kwargs : dict
        The MDP options for the fep calculations on every step. This dictionary must have the structure:

        .. code-block:: text

            {
            'vdw':{
                'step1': <mdp options>,
                'step2': <mdp options>,
                ...
                }
            'coul':{
                'step1': <mdp options>,
                'step2': <mdp options>,
                ...
            'bonded':{
                'step1': <mdp options>,
                'step2': <mdp options>,
                ...
                }
            }
    Raises
    ------
    ValueError
        In case of an invalid ``lambda_type``
    ValueError
        In case of an invalid ``sys_type``
    """
    sim_dir = Path(sim_dir)
    template_dir = Path(template_dir)
    valid_lambda_types = ["vdw", "coul", "bonded"]
    valid_sys_types = ['ligand', 'complex']
    if lambda_type not in valid_lambda_types:
        raise ValueError(f"Non valid lambda_type = {lambda_type}. Must be one of {valid_lambda_types}")
    if sys_type not in valid_sys_types:
        raise ValueError(f"Non valid sys_type = {sys_type}. Must be one of {valid_sys_types}")

    # Take from the source of the package what are the input MDP files
    input_mdp = [step.name for step in list_if_file(template_dir/f"{lambda_type}", ext='.mdp')]

    # Create the lambda string
    lambda_range_str = " ".join(map(str, lambda_values))
    # Create MDP template for fep calculations
    mdp_template = StepMDP(step_path=template_dir/lambda_type)
    for mdp_file in input_mdp:
        step = Path(mdp_file).stem

        # Update MDP step
        mdp_template.set_new_step(step)

        # Check dt and set dt_max if needed, this will be overwrite by the parameters provided in the mdp section of the config
        if 'dt' in mdp_template.parameters:  # Avoid min step, it assumes that the rest of the mdp templates steps have dt defined.
            if float(mdp_template.parameters['dt'].split(';')[0]) > dt_max:
                mdp_template.set_parameters(dt=dt_max)

        # Set, if any, the user MDP options
        if mdp_extra_kwargs:
            try:
                # TODO sanity check on the passed mdp options
                mdp_template.set_parameters(**mdp_extra_kwargs[lambda_type][step])
            except KeyError:
                pass

        # Update lambdas
        mdp_template.set_parameters(**{f"{lambda_type}-lambdas": lambda_range_str})

        # Set to 1 all the bonded-lambdas in case of vdw and coul for the complex
        if sys_type.lower() == 'complex' and lambda_type in ['vdw', 'coul']:
            mdp_template.set_parameters(**{"bonded-lambdas": " ".join(map(str, len(lambda_values)*[1]))})

        for i in range(len(lambda_values)):
            # Create simulation/state/step directory
            (sim_dir/f"simulation/{lambda_type}.{i}/{step}").mkdir(exist_ok=True, parents=True)
            # Update init-lambda-state
            mdp_template.set_parameters(**{"init-lambda-state": i})
            # Write MDP to the proper location
            mdp_template.write(sim_dir/f"simulation/{lambda_type}.{i}/{step}/{step}.mdp")


def get_number_of_frames(input_mdp):
    loaded_mdp_params = MDP().from_file(input_mdp).parameters
    return -(int(loaded_mdp_params['nsteps'].split(';')[0]) // -int(loaded_mdp_params['nstxout-compressed'].split(';')[0])) # ceiling to the next int


if __name__ == "__main__":
    pass
