rule fep_setup_ligand:
    input:
        mdp_vdw=expand(str(TemplatePath.ligand.fep/"vdw/{step}.mdp"), step=[step.stem for step in tools.list_if_file(TemplatePath.ligand.fep/"vdw", ext='.mdp')]),
        mdp_coul=expand(str(TemplatePath.ligand.fep/"coul/{step}.mdp"), step=[step.stem for step in tools.list_if_file(TemplatePath.ligand.fep/"coul", ext='.mdp')])
    params:
        template_dir=str(TemplatePath.ligand.fep),
        vdw_lambdas=config['lambdas']['ligand']['vdw'],
        coul_lambdas=config['lambdas']['ligand']['coul'],
        ligand_names=config['ligand_names'],
        replicas=range(1,1 + config['replicas']),
    output:
        mdp_vdw=expand(out_approach_path+"/{ligand_name}/{replica}/ligand/fep/simulation/vdw.{state}/{step}/{step}.mdp", state=range(len(config['lambdas']['ligand']['vdw'])), step=[step.stem for step in tools.list_if_file(TemplatePath.ligand.fep/"vdw", ext='.mdp')], ligand_name=config['ligand_names'], replica=list(map(str, range(1,1 + config['replicas'])))),
        mdp_coul=expand(out_approach_path+"/{ligand_name}/{replica}/ligand/fep/simulation/coul.{state}/{step}/{step}.mdp", state=range(len(config['lambdas']['ligand']['coul'])), step=[step.stem for step in tools.list_if_file(TemplatePath.ligand.fep/"coul", ext='.mdp')], ligand_name=config['ligand_names'], replica=list(map(str, range(1,1 + config['replicas']))))
    run:

        # In case of user defined MDP keywords, take those from the config
        try:
            # TODO sanity check on the passed MDP options
            mdp_extra_kwargs = config['mdp']['ligand']['fep']
        except KeyError:
            mdp_extra_kwargs = {}

        for ligand_name in params.ligand_names:
            for replica in params.replicas:
                sim_dir = f"{out_approach_path}/{ligand_name}/{replica}/ligand/fep"
                # Create MDP template for Van der Waals states
                mdp.make_fep_dir_structure(
                    sim_dir=sim_dir,
                    template_dir=params.template_dir,
                    lambda_values=params.vdw_lambdas,
                    lambda_type='vdw',
                    sys_type='ligand',
                    dt_max=config["dt_max"],
                    mdp_extra_kwargs=mdp_extra_kwargs,
                )

                # Create MDP template for Coulomb states
                mdp.make_fep_dir_structure(
                    sim_dir=sim_dir,
                    template_dir=params.template_dir,
                    lambda_values=params.coul_lambdas,
                    lambda_type='coul',
                    sys_type='ligand',
                    dt_max=config["dt_max"],
                    mdp_extra_kwargs=mdp_extra_kwargs,
                )