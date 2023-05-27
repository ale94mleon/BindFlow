from abfe.mdp.templates import TemplatePath
from abfe.utils.tools import makedirs, list_if_file
from abfe.mdp import mdp
import os

# Common to all the sub-workflows ligand/replica
run_path = config["run_path"]

rule equil_setup_complex:
    input:
        mdp=expand(TemplatePath.complex.soluble.equi+"/{step}.mdp", step=[os.path.splitext(step)[0] for step in list_if_file(TemplatePath.complex.soluble.equi, ext='mdp')])
    params:
        sim_dir=run_path+"/complex/equil-mdsim",
        template_dir=TemplatePath.complex.soluble.equi,
        # Dynamically access the simulation steps based on the name of the mdp files inside template_dir.
        # Must be defined in this way, outside of the rule is overwrite it.
        steps = [os.path.splitext(step)[0] for step in list_if_file(TemplatePath.complex.soluble.equi, ext='mdp')],
    output:
        mdp=expand(run_path+"/complex/equil-mdsim/{step}/{step}.mdp", step=[os.path.splitext(step)[0] for step in list_if_file(TemplatePath.complex.soluble.equi, ext='mdp')])
    run:
        # Create MDP template
        mdp_template = mdp.StepMDP(step_path = params.template_dir)
        # Not sure if the sorted is needed, just for safety
        for step, input_mdp, output_mdp in sorted(zip(params.steps, input.mdp, output.mdp), key=lambda x: x[0]):
            makedirs(params.sim_dir+f"/{step}")
            # Update MDP step
            mdp_template.set_new_step(step)
            # TODO In case of user defined MDP keywords, take those from the config
            try:
                # TODO sanity check on the passed MDP options
                mdp_template.set_parameters(**config['mdp']['complex']['equi'][step])
            except KeyError:
                pass
            # Write MDP to the proper location
            mdp_template.write(output_mdp)