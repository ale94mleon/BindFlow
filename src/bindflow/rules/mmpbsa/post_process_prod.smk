from bindflow.utils.tools import hard_code_dependencies
from bindflow.utils.tools import gmx_command
from bindflow.utils.tools import load_dependencies_for_shell

import json

approach_path = config["out_approach_path"]
load_dependencies = config['extra_directives']['dependencies']

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


rule create_index_ndx:
    input:
        complex_equi_finished = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.finished",
    output:
        out_index = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/index.ndx",
        out_group = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/groups.json",
    params:
        in_tpr = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.tpr",
    run:
        import logging
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(name)s - %(message)s')
        
        @gmx_command(load_dependencies=hard_code_dependencies + load_dependencies, interactive_script=["q"])
        def make_ndx(**kwargs): ...

        res = make_ndx(f=params.in_tpr, o=output.out_index)
        protein_group = None
        ligand_group = None
        for line in res.stdout.splitlines():
            if protein_group is None:
                if " Protein " in line and " : " in line:
                    protein_group = line.split()[0]
            
            if ligand_group is None:
                if " LIG " in line and " : " in line:
                    ligand_group = line.split()[0]
        
        data = {"protein_group": protein_group, "ligand_group": ligand_group}
        with open(output.out_group, "w") as group_file:
            json.dump(data, group_file)