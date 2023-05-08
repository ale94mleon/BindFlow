from abfe import template
from abfe.utils.tools import makedirs, PathLike
from abfe.utils import mdp
import os
template_dir = template.ligand_fep_template_path

run_path = config["run_path"]
vdw_lambdas = config['ligand_vdw_lambdas']
coul_lambdas = config['ligand_coul_lambdas']

# TODO In case of user defined MDP keywords, take those from the config
try:
    # TODO sanity check on the passed MDP options
    mdp_extra_kwargs = config['mdp']['ligand']['fep']
except KeyError:
    mdp_extra_kwargs = {}

vdw_steps = [os.path.splitext(step)[0] for step in list_if_file(template_dir+"/vdw", ext='mdp')]
coul_steps = [os.path.splitext(step)[0] for step in list_if_file(template_dir+"/coul", ext='mdp')]

rule fep_setup_ligand:
    input:
        mdp_vdw=expand(template_dir+"/vdw/{step}.mdp", step=vdw_steps),
        mdp_coul=expand(template_dir+"/coul/{step}.mdp", step=coul_steps)
    params:
        sim_dir=run_path+"/ligand/fep",
        vdw_range_str=" ".join(map(str, vdw_lambdas)),
        coul_range_str=" ".join(map(str, coul_lambdas)),
    output:
        mdp_vdw=expand(run_path+"/ligand/fep/simulation/vdw.{state}/{step}/{step}.mdp", state=range(len(vdw_lambdas)), step=vdw_steps),
        mdp_coul=expand(run_path+"/ligand/fep/simulation/coul.{state}/{step}/{step}.mdp", state=range(len(coul_lambdas)), step=coul_steps)
    run:
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