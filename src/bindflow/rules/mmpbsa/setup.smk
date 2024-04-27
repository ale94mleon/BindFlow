import bindflow.mmpbsa_in.mmpbsa_in_load_defaults

approach_path = config["out_approach_path"]
load_dependencies = config['extra_directives']['dependencies']


# TODO move to input preparation at the beginning, this go with the topology, configuration and index files
# The parameters should be expose to user customization. An Class like the MDP class may do the job
# Then the new parameters are read from the yml file and update the default values, forcefields and other should be protected and not mutable for MMPBSA (force fields are 
# generated at the beginning.)
rule create_mmpbsa_in:
    input:
    output:
        mmpbsa_in_file = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/mmpbsa.in",
    run:
        import logging
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(name)s - %(message)s')
        
        mmpbsa_config = {}
        if "mmpbsa" in config.keys():
            mmpbsa_config = config["mmpbsa"]
        if not "pb" in mmpbsa_config.keys() and not "gb" in mmpbsa_config.keys():
            mmpbsa_config["pb"] = True

        input_file = bindflow.mmpbsa_in.mmpbsa_in_load_defaults.gmxmmpbsa_input_file(mmpbsa_in_opts = mmpbsa_config)
        input_file.store_template(output.mmpbsa_in_file)