from abfe import template
from abfe.utils.tools import makedirs, list_if_file
from abfe.utils import mdp
import os, shutil

# Common to all the sub-workflows ligand/replica
input_path = config['input_data_path']
run_path = config["run_path"]

rule fep_setup_complex:
    input:
        complex_top=input_path+"/complex/complex.top",
        boresch_top=run_path+"/complex/equil-mdsim/boreschcalc/BoreschRestraint.top",
        mdp_vdw=expand(template.complex_fep_template_path+"/vdw/{step}.mdp", step=[os.path.splitext(step)[0] for step in list_if_file(template.complex_fep_template_path+"/vdw", ext='mdp')]),
        mdp_coul=expand(template.complex_fep_template_path+"/coul/{step}.mdp", step=[os.path.splitext(step)[0] for step in list_if_file(template.complex_fep_template_path+"/coul", ext='mdp')]),
        mdp_bonded=expand(template.complex_fep_template_path+"/bonded/{step}.mdp", step=[os.path.splitext(step)[0] for step in list_if_file(template.complex_fep_template_path+"/bonded", ext='mdp')])
    params:
        sim_dir=run_path+"/complex/fep",
        template_dir=template.complex_fep_template_path,
        vdw_lambdas=config['lambdas']['complex']['vdw'],
        coul_lambdas=config['lambdas']['complex']['coul'],
        bonded_lambdas=config['lambdas']['complex']['bonded'],
    output:
        mdp_vdw=expand(run_path+"/complex/fep/simulation/vdw.{state}/{step}/{step}.mdp", state=range(len(config['lambdas']['complex']['vdw'])), step=[os.path.splitext(step)[0] for step in list_if_file(template.complex_fep_template_path+"/vdw", ext='mdp')]),
        mdp_coul=expand(run_path+"/complex/fep/simulation/coul.{state}/{step}/{step}.mdp", state=range(len(config['lambdas']['complex']['coul'])), step=[os.path.splitext(step)[0] for step in list_if_file(template.complex_fep_template_path+"/coul", ext='mdp')]),
        mdp_bonded=expand(run_path+"/complex/fep/simulation/bonded.{state}/{step}/{step}.mdp", state=range(len(config['lambdas']['complex']['bonded'])), step=[os.path.splitext(step)[0] for step in list_if_file(template.complex_fep_template_path+"/bonded", ext='mdp')]),
        fep_top=run_path+"/complex/fep/fep-topology/complex_boresch.top",
    run:

        # TODO In case of user defined MDP keywords, take those from the config
        try:
            # TODO sanity check on the passed MDP options
            mdp_extra_kwargs = config['mdp']['complex']['fep']
        except KeyError:
            mdp_extra_kwargs = {}
        
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
            template_dir = params.template_dir,
            lambda_values = params.vdw_lambdas,
            lambda_type = 'vdw',
            **mdp_extra_kwargs,
        )

        # Create MDP template for Coulomb states
        mdp.make_fep_dir_structure(
            sim_dir = params.sim_dir,
            template_dir = params.template_dir,
            lambda_values = params.coul_lambdas,
            lambda_type = 'coul',
            **mdp_extra_kwargs,
        )

        # Create MDP template for Restraint states
        mdp.make_fep_dir_structure(
            sim_dir = params.sim_dir,
            template_dir = params.template_dir,
            lambda_values = params.bonded_lambdas,
            lambda_type = 'bonded',
            **mdp_extra_kwargs,
        )