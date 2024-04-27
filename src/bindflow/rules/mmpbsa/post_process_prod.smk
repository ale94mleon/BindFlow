from bindflow.utils.tools import hard_code_dependencies
from bindflow.utils.tools import gmx_command
from bindflow.utils.tools import load_dependencies_for_shell

import json

approach_path = config["out_approach_path"]
load_dependencies = config['extra_directives']['dependencies']


# TODO: use the same scheme used for boresch restraint to fix the trajectory, maybe inside the rule of gmx_run, 
rule remove_PBC_in_trajectory:
    input:
        complex_equi_finished = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.finished",
    output:
        out_xtc = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod_noPBC.xtc",
    params:
        in_tpr = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.tpr",
        in_xtc = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.xtc",
    run:
        import logging
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(name)s - %(message)s')
        
        @gmx_command(load_dependencies=hard_code_dependencies + load_dependencies, interactive_script=["1", "0"])
        def trjconv(**kwargs): ...

        trjconv(s=params.in_tpr, f=params.in_xtc, o=output.out_xtc, pbc="mol", center=True, ur="compact")