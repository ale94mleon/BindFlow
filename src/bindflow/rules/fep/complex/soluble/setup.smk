rule fep_setup_complex:
    input:
        complex_top_dir=out_approach_path+"/{ligand_name}/input/complex",
        boresch_top=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/boreschcalc/BoreschRestraint.top",
        mdp_vdw=expand(str(TemplatePath.complex.soluble.fep/"vdw/{step}.mdp"), step=[step.stem for step in tools.list_if_file(TemplatePath.complex.soluble.fep/"vdw", ext='.mdp')]),
        mdp_coul=expand(str(TemplatePath.complex.soluble.fep/"coul/{step}.mdp"), step=[step.stem for step in tools.list_if_file(TemplatePath.complex.soluble.fep/"coul", ext='.mdp')]),
        mdp_bonded=expand(str(TemplatePath.complex.soluble.fep/"bonded/{step}.mdp"), step=[step.stem for step in tools.list_if_file(TemplatePath.complex.soluble.fep/"bonded", ext='.mdp')])
    params:
        template_dir=TemplatePath.complex.soluble.fep,
        vdw_lambdas=config['lambdas']['complex']['vdw'],
        coul_lambdas=config['lambdas']['complex']['coul'],
        bonded_lambdas=config['lambdas']['complex']['bonded'],
    output:
        mdp_vdw=expand(out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/vdw.{state}/{step}/{step}.mdp", state=range(len(config['lambdas']['complex']['vdw'])), step=[step.stem for step in tools.list_if_file(TemplatePath.complex.soluble.fep/"vdw", ext='.mdp')], allow_missing=True),
        mdp_coul=expand(out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/coul.{state}/{step}/{step}.mdp", state=range(len(config['lambdas']['complex']['coul'])), step=[step.stem for step in tools.list_if_file(TemplatePath.complex.soluble.fep/"coul", ext='.mdp')], allow_missing=True),
        mdp_bonded=expand(out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/bonded.{state}/{step}/{step}.mdp", state=range(len(config['lambdas']['complex']['bonded'])), step=[step.stem for step in tools.list_if_file(TemplatePath.complex.soluble.fep/"bonded", ext='.mdp')], allow_missing=True),
        fep_top=out_approach_path+"/{ligand_name}/{replica}/complex/fep/topology/complex_boresch.top",
    run:
        # In case of user defined MDP keywords, take those from the config
        try:
            # TODO sanity check on the passed MDP options
            mdp_extra_kwargs = config['mdp']['complex']['fep']
        except KeyError:
            mdp_extra_kwargs = {}
        
        sim_dir = Path(out_approach_path)/f"{wildcards.ligand_name}/{wildcards.replica}/complex/fep"
        
        # Make topology directory
        out_top_dir = sim_dir/"topology"
        out_top_dir.mkdir(exist_ok=True, parents=True)
        
        # Copy complex topology (only itp),
        # I can not set them as output of the rule becasue the name might change
        for itp_file in tools.list_if_file(input.complex_top_dir, ext='.itp'):
            shutil.copy(itp_file, out_top_dir)
        
        # Modify the main topology incorporating the boresch restraints
        with open(Path(input.complex_top_dir)/'complex.top', 'r') as original_top:
            with open(input.boresch_top, 'r') as boresch_top:
                with open(output.fep_top, 'w') as final_top:
                    final_top.write(original_top.read() + boresch_top.read())

        # Create MDP template for Van der Waals states
        mdp.make_fep_dir_structure(
            sim_dir=sim_dir,
            template_dir=params.template_dir,
            lambda_values=params.vdw_lambdas,
            lambda_type='vdw',
            sys_type='complex',
            dt_max=config["dt_max"],
            mdp_extra_kwargs=mdp_extra_kwargs,
        )

        # Create MDP template for Coulomb states
        mdp.make_fep_dir_structure(
            sim_dir=sim_dir,
            template_dir=params.template_dir,
            lambda_values=params.coul_lambdas,
            lambda_type='coul',
            sys_type='complex',
            dt_max=config["dt_max"],
            mdp_extra_kwargs=mdp_extra_kwargs,
        )

        # Create MDP template for Restraint states
        mdp.make_fep_dir_structure(
            sim_dir=sim_dir,
            template_dir=params.template_dir,
            lambda_values=params.bonded_lambdas,
            lambda_type='bonded',
            sys_type='complex',
            dt_max=config["dt_max"],
            mdp_extra_kwargs=mdp_extra_kwargs,
        )