import bindflow.mmpbsa_in.mmpbsa_in_load_defaults
from bindflow.mdp.templates import mdp
from bindflow.mdp.templates import TemplatePath
from bindflow.utils import tools
import os
approach_path = config["out_approach_path"]


# TODO
# The parameters should be expose to user customization. A Class like the MDP class may do the job
# Then the new parameters are read from the yml file and update the default values, forcefields and other should be protected and not mutable for MMPBSA (force fields are 
# generated at the beginning.)
rule create_mmpbsa_in:
    output:
        mmpbsa_in = approach_path + "/{ligand_name}/input/mmpbsa.in",
    run:
        import logging
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(name)s - %(message)s')
        
        mmpbsa_config = {}
        if "mmpbsa" in config.keys():
            mmpbsa_config = config["mmpbsa"]
        if not "pb" in mmpbsa_config.keys() and not "gb" in mmpbsa_config.keys():
            mmpbsa_config["pb"] = True

        input_file = bindflow.mmpbsa_in.mmpbsa_in_load_defaults.gmxmmpbsa_input_file(mmpbsa_in_opts = mmpbsa_config)
        input_file.store_template(output.mmpbsa_in)

rule mmpbsa_setup:
    output:
        mdp = expand(approach_path + "/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{split}/prod.mdp", split = range(config['split']), ligand_name = config['ligand_names'], replica = list(map(str, range(1,1 + config['replicas']))))
    run:
        if config["complex_type"] == 'soluble':
            template_dir = TemplatePath.complex.soluble.mmpbsa
        elif config["complex_type"] =='membrane':
            template_dir = TemplatePath.complex.membrane.mmpbsa
        mdp_template = mdp.MDP(template_dir / 'prod.mdp').from_file()
        # In case of user defined MDP keywords, take those from the config
        try:
            # TODO sanity check on the passed MDP options
            mdp_template.set_parameters(**config['mdp']['complex']['mmpbsa']['prod'])
        except KeyError:
            pass
        
        for m in output.mdp:
            tools.makedirs(os.path.driname(m))
            mdp_template.write(m)
        