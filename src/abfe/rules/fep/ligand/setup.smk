from abfe.utils.tools import makedirs, list_if_file
from abfe.mdp.templates import TemplatePath
from abfe.mdp import mdp
import os


# Common to all the sub-workflows ligand/replica
run_path = config["run_path"]


rule fep_setup_ligand:
    input:
        mdp_vdw=expand(TemplatePath.ligand.fep+"/vdw/{step}.mdp", step=[os.path.splitext(step)[0] for step in list_if_file(TemplatePath.ligand.fep+"/vdw", ext='mdp')]),
        mdp_coul=expand(TemplatePath.ligand.fep+"/coul/{step}.mdp", step=[os.path.splitext(step)[0] for step in list_if_file(TemplatePath.ligand.fep+"/coul", ext='mdp')])
    params:
        sim_dir=run_path+"/ligand/fep",
        template_dir = TemplatePath.ligand.fep,
        vdw_lambdas = config['lambdas']['ligand']['vdw'],
        coul_lambdas = config['lambdas']['ligand']['coul'],

    output:
        mdp_vdw=expand(run_path+"/ligand/fep/simulation/vdw.{state}/{step}/{step}.mdp", state=range(len(config['lambdas']['ligand']['vdw'])), step=[os.path.splitext(step)[0] for step in list_if_file(TemplatePath.ligand.fep+"/vdw", ext='mdp')]),
        mdp_coul=expand(run_path+"/ligand/fep/simulation/coul.{state}/{step}/{step}.mdp", state=range(len(config['lambdas']['ligand']['coul'])), step=[os.path.splitext(step)[0] for step in list_if_file(TemplatePath.ligand.fep+"/coul", ext='mdp')])
    run:

        # TODO In case of user defined MDP keywords, take those from the config
        try:
            # TODO sanity check on the passed MDP options
            mdp_extra_kwargs = config['mdp']['ligand']['fep']
        except KeyError:
            mdp_extra_kwargs = {}

        # Create MDP template for Van der Waals states
        mdp.make_fep_dir_structure(
            sim_dir = params.sim_dir,
            template_dir = params.template_dir,
            lambda_values = params.vdw_lambdas,
            lambda_type = 'vdw',
            sys_type = 'ligand',
            mdp_extra_kwargs = mdp_extra_kwargs,
        )

        # Create MDP template for Coulomb states
        mdp.make_fep_dir_structure(
            sim_dir = params.sim_dir,
            template_dir = params.template_dir,
            lambda_values = params.coul_lambdas,
            lambda_type = 'coul',
            sys_type = 'ligand',
            mdp_extra_kwargs = mdp_extra_kwargs,
        )