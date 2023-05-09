from abfe import template
from abfe.utils.tools import makedirs, list_if_file
from abfe.utils import mdp
import os, shutil
template_dir=template.complex_fep_template_path

input_path = config['input_data_path']
run_path = config["run_path"]
coul_lambdas = config['lambdas']['complex']['coul']
vdw_lambdas = config['lambdas']['complex']['vdw']
bonded_lambdas = config['lambdas']['complex']['bonded']

vdw_steps = [os.path.splitext(step)[0] for step in list_if_file(template_dir+"/vdw", ext='mdp')]
coul_steps = [os.path.splitext(step)[0] for step in list_if_file(template_dir+"/coul", ext='mdp')]
bonded_steps = [os.path.splitext(step)[0] for step in list_if_file(template_dir+"/bonded", ext='mdp')]

# TODO In case of user defined MDP keywords, take those from the config
try:
    # TODO sanity check on the passed MDP options
    mdp_extra_kwargs = config['mdp']['complex']['fep']
except KeyError:
    mdp_extra_kwargs = {}

rule fep_setup_complex:
    input:
        complex_top=input_path+"/complex/complex.top",
        boresch_top=run_path+"/complex/equil-mdsim/boreschcalc/BoreschRestraint.top",
        mdp_vdw=expand(template_dir+"/vdw/{step}.mdp", step=vdw_steps),
        mdp_coul=expand(template_dir+"/coul/{step}.mdp", step=coul_steps),
        mdp_bonded=expand(template_dir+"/bonded/{step}.mdp", step=bonded_steps)
    params:
        sim_dir=run_path+"/complex/fep",
    output:
        mdp_vdw=expand(run_path+"/complex/fep/simulation/vdw.{state}/{step}/{step}.mdp", state=range(len(vdw_lambdas)), step=vdw_steps),
        mdp_coul=expand(run_path+"/complex/fep/simulation/coul.{state}/{step}/{step}.mdp", state=range(len(coul_lambdas)), step=coul_steps),
        mdp_bonded=expand(run_path+"/complex/fep/simulation/bonded.{state}/{step}/{step}.mdp", state=range(len(bonded_lambdas)), step=coul_steps),
        fep_top=run_path+"/complex/fep/fep-topology/complex_boresch.top",
    run:
        # Copy complex topology
        shutil.copytree(input.complex_top, sim_dir+"/complex/fep/fep-topology")
        
        # Modify the main topology incorporating the boresch restraints
        with open(params.sim_dir+"/fep-topology/complex.top", 'r') as original_top:
            with open(input.boresch_top, 'r') as boresch_top:
                with open(output.fep_top, 'w') as final_top:
                    final_top.write(original_top.read() + boresch_top.read())

        # Create MDP template for Van der Waals states
        mdp.make_fep_dir_structure(
            sim_dir = params.sim_dir,
            template_dir = template_dir,
            lambda_values = vdw_lambdas,
            lambda_type = 'vdw',
            **mdp_extra_kwargs,
        )

        # Create MDP template for Coulomb states
        mdp.make_fep_dir_structure(
            sim_dir = params.sim_dir,
            template_dir = template_dir,
            lambda_values = coul_lambdas,
            lambda_type = 'coul',
            **mdp_extra_kwargs,
        )

        # Create MDP template for Restraint states
        mdp.make_fep_dir_structure(
            sim_dir = params.sim_dir,
            template_dir = template_dir,
            lambda_values = coul_lambdas,
            lambda_type = 'bonded',
            **mdp_extra_kwargs,
        )