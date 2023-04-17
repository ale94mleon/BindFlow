import json, os
from typing import Union
from abfe import template as tp

PathLike = Union[os.PathLike, str, bytes]

_MDP_MEMB_EQUI_TEMPLATES_PATH = os.path.join(tp.membrane_path,'complex_equil_workflow')
_MDP_PARAM_DEFAULT={
    "integrator": "steep",
    "emtol": "1000.0",
    "nsteps": "5000",
    "nstlist": "10",
    "cutoff-scheme": "Verlet",
    "rlist": "1.0",
    "vdwtype": "Cut-off",
    "vdw-modifier": "Potential-shift-Verlet",
    "rvdw_switch": "0",
    "rvdw": "1.0",
    "coulombtype": "pme",
    "rcoulomb": "1.0",
    "epsilon-r": "1",
    "epsilon-rf": "1",
    "constraints": "h-bonds",
    "constraint_algorithm": "LINCS"
}
class MDP:
    def __init__(self, **kwargs):
        self.parameters = {}
        self._set_default_parameters()
        self.set_parameters(**kwargs)

    def _set_default_parameters(self):
        # Add any default parameters here
        self.parameters = _MDP_PARAM_DEFAULT

    def set_parameters(self, **kwargs):
        self.parameters.update(kwargs)

    def from_file(self, template_filename, clean_current_parameters = True):
        with open(template_filename, 'r') as f:
            lines = f.readlines()
            if clean_current_parameters:
                self.parameters = {}
            for line in lines:
                if line.startswith(';') or line.startswith('#'):
                    continue
                tokens = line.strip().split('=', 1)
                if len(tokens) != 2:
                    continue
                parameter_name = tokens[0].strip()
                parameter_value = tokens[1].strip()
                self.parameters[parameter_name] = parameter_value
    
    def to_string(self):
        s = ''
        for parameter_name, parameter_value in self.parameters.items():
            s += f'{parameter_name:<40} = {parameter_value}\n'
        return s

    def write(self, filename:str):
        with open(filename, 'w') as f:
            f.write(self.to_string())

    def __repr__(self):
        return f"{self.__class__.__name__}({json.dumps(self.parameters, indent=4)})"

class MembraneEquiMDP(MDP):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
    def from_archive(self, name = 'step6.0_minimization'):
        valid_names = ['step6.0_minimization'] + [f'step6.{i}_equilibration' for i in range(1,8)]
        if name not in valid_names:
            raise ValueError(f"name = {name} is not valid, must be one of: {valid_names}")
        self.from_file(os.path.join(_MDP_MEMB_EQUI_TEMPLATES_PATH, f"{name}.mdp"))


eq_mdp = MembraneEquiMDP()
eq_mdp.from_archive(name = 'step6.7_equilibration') # name = 'step6.5_equilibration'
print(eq_mdp)