from bindflow.utils import tools
from bindflow.preparation import boresch
from bindflow.mdp import mdp


rule equil_complex_get_boresch_restraints:
    input:
        finished = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.finished",
        mdp = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.mdp",
    params:
        in_tpr = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.tpr",
        in_xtc = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.xtc",
        run_dir = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/boreschcalc",
    output:
        gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/boreschcalc/ClosestRestraintFrame.gro",
        top = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/boreschcalc/BoreschRestraint.top",
        boresch_dG_off = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/boreschcalc/dG_off.dat"
    run:
        # Set default host name (reachable as gmx index) and selection (valid MDAnalysis selection)
        host_name = 'Protein'
        host_selection = 'protein and name CA'
        # Only if debug is activated, update name and selection based on environment variable
        if 'abfe_debug' in os.environ:
            if os.environ['abfe_debug'] == 'True': # environ save the variables as strings 
                if 'abfe_debug_host_name' in os.environ:
                    host_name = os.environ['abfe_debug_host_name']
                if 'abfe_debug_host_selection' in os.environ:
                    host_selection = os.environ['abfe_debug_host_selection']
        # Fix trajectory.
        tools.makedirs(params.run_dir)
        tools.run(f"export GMX_MAXBACKUP=-1; echo 'System' | gmx trjconv -s {params.in_tpr} -f {params.in_xtc} -o {params.run_dir}/whole.xtc -pbc whole")
        tools.run(f"export GMX_MAXBACKUP=-1; echo 'System' | gmx trjconv -s {params.in_tpr} -f {params.run_dir}/whole.xtc -o {params.run_dir}/nojump.xtc -pbc nojump")
        tools.run(f"export GMX_MAXBACKUP=-1; echo '{host_name} System' | gmx trjconv -s {params.in_tpr} -f {params.run_dir}/nojump.xtc -o {params.run_dir}/prod_center.xtc -pbc mol -center -ur compact")
        # Clean
        tools.run(f"rm {params.run_dir}/whole.xtc {params.run_dir}/nojump.xtc")

        # Getting Borech restraints
        mdp_params = mdp.MDP().from_file(input.mdp).parameters
        if 'ref-t' in mdp_params:
            temperature = float(mdp_params['ref-t'].split()[0])
        elif 'ref_t' in mdp_params:
            temperature = float(mdp_params['ref_t'].split()[0])
        boresch.gen_restraint(
            topology = params.in_tpr,
            trajectory = f"{params.run_dir}/prod_center.xtc",
            outpath = params.run_dir,
            temperature = temperature,
            host_selection = host_selection
        )
        # Clean
        tools.run(f"rm {params.run_dir}/prod_center.xtc")
