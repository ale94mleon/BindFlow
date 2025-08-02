rule equil_setup_complex:
    input:
        mdp=expand(str(TemplatePath.complex.membrane.equi/"{step}.mdp"), step=[step.stem for step in tools.list_if_file(TemplatePath.complex.membrane.equi, ext='.mdp')])
    params:
        template_dir=TemplatePath.complex.membrane.equi,
        # Dynamically access the simulation steps based on the name of the mdp files inside template_dir.
        # Must be defined in this way, outside of the rule is overwrite it.
        steps=[step.stem for step in tools.list_if_file(TemplatePath.complex.membrane.equi, ext='.mdp')],
        ligand_names=config['ligand_names'],
        replicas=range(1,1 + config['replicas']),
    output:
        mdp=expand(out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/{step}/{step}.mdp", step=[step.stem for step in tools.list_if_file(TemplatePath.complex.membrane.equi, ext='.mdp')], ligand_name=config['ligand_names'], replica=list(map(str, range(1,1 + config['replicas']))))
    run:
       for ligand_name in params.ligand_names:
            for replica in params.replicas:
                sim_dir=Path(out_approach_path)/f"{ligand_name}/{replica}/complex/equil-mdsim"
                # Create MDP template

                mdp_template=mdp.StepMDP(step_path = params.template_dir)

                # Not sure if the sorted is needed, just for safety
                for step in params.steps:
                    (sim_dir/step).mkdir(exist_ok=True, parents=True)
                    output_mdp = sim_dir/step/f"{step}.mdp"
                    
                    # Change to new MDP step
                    mdp_template.set_new_step(step)
                    
                    # Check dt and set dt_max if needed, this will be overwrite by the parameters provided in the mdp section of the config
                    if 'dt' in mdp_template.parameters: # Avoid min step, it assumes that the rest of the mdp templates steps have dt defined.
                        if float(mdp_template.parameters['dt'].split(';')[0]) > config["dt_max"]:
                            mdp_template.set_parameters(dt=config["dt_max"])
                        
                    # In case of user defined MDP keywords, take those from the config
                    try:
                        # TODO sanity check on the passed MDP options
                        mdp_template.set_parameters(**config['mdp']['complex']['equi'][step])
                    except KeyError:
                        pass
                    # Write MDP to the proper location
                    mdp_template.write(output_mdp)