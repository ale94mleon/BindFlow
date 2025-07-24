from bindflow.mmpbsa_in.input_loader import GMXMMPBSAInputMaker
import tempfile

# Common to all sub-workflows
ligand_names = config["ligand_names"]
samples = list(map(str, range(1,1 + config["samples"])))
replicas = list(map(str, range(1,1 + config["replicas"])))

rule mmxbsa_setup:
    input:
        finished=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.finished",
        mdp=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.mdp"
    output:
        mdp=expand(out_approach_path+"/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/prod.mdp", sample=samples, allow_missing=True),
        gro=expand(out_approach_path+"/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/init.gro", sample=samples, allow_missing=True)
    params:
        in_tpr=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.tpr",
        in_xtc=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.xtc",
        sim_dir=out_approach_path+"/{ligand_name}/{replica}/complex/mmpbsa/simulation",
    run:
        # Update MDP with user provided options
        if config["complex_type"] == 'soluble':
            template_dir = TemplatePath.complex.soluble.mmpbsa
        elif config["complex_type"] == 'membrane':
            template_dir = TemplatePath.complex.membrane.mmpbsa
        mdp_template = mdp.MDP().from_file(template_dir/'prod.mdp')
        # In case of user defined MDP keywords, take those from the config
        try:
            # TODO sanity check on the passed MDP options
            mdp_template.set_parameters(**config['mdp']['complex']['mmpbsa']['prod'])
        except KeyError:
            pass
        
        sim_dir = Path(params.sim_dir).resolve()
        sim_dir.mkdir(exist_ok=True, parents=True)

        skip = 1

        with tempfile.TemporaryDirectory(prefix='split_', dir=sim_dir) as tmp_dir:
            @tools.gmx_command(load_dependencies=load_dependencies, stdin_command="echo \"System\"")
            def trjconv(**kwargs): ...

            if skip > 0:
                trjconv(f=params.in_xtc, s=params.in_tpr, o=f"{tmp_dir}/.gro", sep=True, skip=skip)
            else:
                trjconv(f=params.in_xtc, s=params.in_tpr, o=f"{tmp_dir}/.gro", sep=True)

            frames = list(Path(tmp_dir).glob('*.gro'))
            if len(frames) < len(output.gro):
                raise RuntimeError("Not enough frames in equil-mdsim/prod/prod.xtc")
            frames = tools.natsort(frames)

            for frame, gro, mdp_file in zip(frames, output.gro, output.mdp):
                Path(gro).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy(frame, gro)
                mdp_template.write(mdp_file)
            print(f"Generated a total of {frames} frames. Using the first {len(output.gro)} frames for MM(P/G)BSA production simulations.")


rule create_mmxbsa_in:
    input:
        mdp=expand(out_approach_path+"/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/prod.mdp", sample=samples, replica=replicas, allow_missing=True)
    output:
        mmpbsa_in=out_approach_path+"/{ligand_name}/input/mmpbsa.in"
    run:
        if "mmpbsa" in config.keys():
            mmpbsa_config = config["mmpbsa"]
        else:
            mmpbsa_config = {}
        
        # Get simulation temperature from any MDP file (all should have the same)
        mdp_params = mdp.MDP().from_file(input.mdp[0]).parameters
        if 'ref-t' in mdp_params:
            temperature = float(mdp_params['ref-t'].split()[0])
        elif 'ref_t' in mdp_params:
            temperature = float(mdp_params['ref_t'].split()[0])
        
        if 'general' in mmpbsa_config:
            mmpbsa_config['general']['temperature'] = temperature
        else:
            mmpbsa_config['general'] = {"temperature":temperature}
        
        # FIXME:
        # Ion concentration is used for gb calculation, this is set during solvation step, not sure if needed
        # There is also a membrane parameter 
        
        input_file = GMXMMPBSAInputMaker(**mmpbsa_config)
        
        input_file.write(output.mmpbsa_in)