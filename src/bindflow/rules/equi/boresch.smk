from bindflow.preparation import boresch

rule equil_complex_get_boresch_restraints:
    input:
        finished=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.finished",
        mdp=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.mdp",
    params:
        in_tpr=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.tpr",
        in_xtc=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.xtc",
        run_dir=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/boreschcalc",
    output:
        gro=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/boreschcalc/ClosestRestraintFrame.gro",
        top=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/boreschcalc/BoreschRestraint.top",
        boresch_dG_off=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/boreschcalc/dG_off.dat"
    run:
        # Set default host name (reachable as gmx index) and selection (valid MDAnalysis selection)
        host_name='Protein'
        host_selection='protein and name CA'
        # Only if debug is activated, update name and selection based on environment variable
        if 'abfe_debug' in os.environ:
            if os.environ['abfe_debug'] == 'True': # environ save the variables as strings 
                if 'abfe_debug_host_name' in os.environ:
                    host_name = os.environ['abfe_debug_host_name']
                if 'abfe_debug_host_selection' in os.environ:
                    host_selection = os.environ['abfe_debug_host_selection']

        # Fix trajectory.
        tools.center_xtc(
            tpr=params.in_tpr,
            xtc=params.in_xtc,
            run_dir=params.run_dir,
            host_name=host_name
        )

        # Getting Borech restraints
        mdp_params = mdp.MDP().from_file(input.mdp).parameters
        if 'ref-t' in mdp_params:
            temperature = float(mdp_params['ref-t'].split()[0])
        elif 'ref_t' in mdp_params:
            temperature = float(mdp_params['ref_t'].split()[0])
        boresch.gen_restraint(
            topology=params.in_tpr,
            trajectory=f"{params.run_dir}/center.xtc",
            outpath=params.run_dir,
            temperature=temperature,
            host_selection=host_selection
        )
        # Clean
        (Path(params.run_dir)/"center.xtc").unlink()