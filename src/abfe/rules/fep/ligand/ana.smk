from abfe.scripts.free_energy import analysis
from abfe.utils import tools
from abfe.mdp import mdp
import os

run_path = config["run_path"]
threads = config['threads']
# Ana
rule fep_ana_get_dg_ligand_contributions:
    input:
        xvg_vdw_loc=expand(run_path+"/ligand/fep/simulation/vdw.{state}/prod/prod.xvg", state=range(len(config['lambdas']['ligand']['vdw']))),
        xvg_coul_loc=expand(run_path+"/ligand/fep/simulation/coul.{state}/prod/prod.xvg", state=range(len(config['lambdas']['ligand']['coul']))),
        # Make sure that the simualtion ends properly
        gro_vdw_loc=expand(run_path+"/ligand/fep/simulation/vdw.{state}/prod/prod.gro", state=range(len(config['lambdas']['ligand']['vdw']))),
        gro_coul_loc=expand(run_path+"/ligand/fep/simulation/coul.{state}/prod/prod.gro", state=range(len(config['lambdas']['ligand']['coul']))),
        # To get the simulaiton temperature
        mdp_vdw_0_prod=run_path+"/ligand/fep/simulation/vdw.0/prod/prod.mdp",
        # TODO: This dependency makes wait the ligand simulations. Could be used on the complex/ana???
        boresch_dat = run_path+"/complex/equil-mdsim/boreschcalc/dG_off.dat",
    params:
        ana_loc=run_path+"/ligand/fep/ana",
    output:
        ligand_json=run_path+"/ligand/fep/ana/dg_ligand_contributions.json"
    threads: threads # TODO: Sometimes the rule hang for a long time
    run:
        # Make directory
        tools.makedirs(params.ana_loc)        
        # Get the simulaiton temperature from the prod.mdp of the state 0 of vdw
        mdp_params = mdp.MDP().from_file(input.mdp_vdw_0_prod).parameters
        if 'ref-t' in mdp_params:
            temperature = float(mdp_params['ref-t'].split()[0])
        elif 'ref_t' in mdp_params:
            temperature = float(mdp_params['ref_t'].split()[0])
        
        analysis.get_dG_contributions(
            boresch_data = input.boresch_dat,
            out_json_path = output.ligand_json,
            # Check if it is necessary to remove some initial burning simulation time
            lower = None,
            upper = None,
            min_samples = 500,
            temperature = temperature,
            # Sort the paths
            vdw = sorted(input.xvg_vdw_loc, key=lambda x: int(os.path.normpath(x).split(os.path.sep)[-3].split('.')[-1])),
            coul = sorted(input.xvg_coul_loc, key=lambda x: int(os.path.normpath(x).split(os.path.sep)[-3].split('.')[-1])),
        )
