from abfe.scripts.free_energy import analysis
from abfe.utils import tools
from abfe.mdp import mdp
import os
import shutil

run_path = config["run_path"]
threads = config['threads']

# Ana
rule fep_ana_get_dg_complex_contributions:
    input:
        xvg_vdw_loc=expand(run_path+"/complex/fep/simulation/vdw.{state}/prod/prod.xvg", state=range(len(config['lambdas']['complex']['vdw']))),
        xvg_coul_loc=expand(run_path+"/complex/fep/simulation/coul.{state}/prod/prod.xvg", state=range(len(config['lambdas']['complex']['coul']))),
        xvg_bonded_loc=expand(run_path+"/complex/fep/simulation/bonded.{state}/prod/prod.xvg", state=range(len(config['lambdas']['complex']['bonded']))),
        # Make sure that the simualtion ends properly
        gro_vdw_loc=expand(run_path+"/complex/fep/simulation/vdw.{state}/prod/prod.gro", state=range(len(config['lambdas']['complex']['vdw']))),
        gro_coul_loc=expand(run_path+"/complex/fep/simulation/coul.{state}/prod/prod.gro", state=range(len(config['lambdas']['complex']['coul']))),
        gro_bonded_loc=expand(run_path+"/complex/fep/simulation/bonded.{state}/prod/prod.gro", state=range(len(config['lambdas']['complex']['bonded']))),
        # To get the simulaiton temperature
        mdp_vdw_0_prod=run_path+"/complex/fep/simulation/vdw.0/prod/prod.mdp",
    params:
        ana_loc=run_path+"/complex/fep/ana",
    output:
        complex_json=run_path+"/complex/fep/ana/dg_complex_contributions.json"
    threads: threads # TODO: Sometimes the rule hang for a long time
    run:
        # Make directory
        tools.makedirs(params.ana_loc)
        # Get the simulaiton temperature from the prod.mdp of the state 0 of vdw
        temperature = float(mdp.MDP().from_file(input.mdp_vdw_0_prod).parameters['ref-t'].split()[0]) 
        analysis.get_dG_contributions(
            boresch_data = None,
            out_json_path = output.complex_json,
            # Check if it is necessary to remove some initial burning simulation time
            lower = None,
            upper = None,
            min_samples = 500,
            temperature = temperature,
            # Sort the paths
            vdw = sorted(input.xvg_vdw_loc, key=lambda x: int(os.path.normpath(x).split(os.path.sep)[-3].split('.')[-1])),
            coul = sorted(input.xvg_coul_loc, key=lambda x: int(os.path.normpath(x).split(os.path.sep)[-3].split('.')[-1])),
            bonded = sorted(input.xvg_bonded_loc, key=lambda x: int(os.path.normpath(x).split(os.path.sep)[-3].split('.')[-1])),
        )

