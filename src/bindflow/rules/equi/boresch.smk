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
        # Fix trajectory.
        tools.center_xtc(
            tpr=params.in_tpr,
            xtc=params.in_xtc,
            run_dir=params.run_dir,
            host_name=config["host_name"],
            load_dependencies=load_dependencies
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
            host_selection=config["host_selection"]
        )
        # Clean
        (Path(params.run_dir)/"center.xtc").unlink()