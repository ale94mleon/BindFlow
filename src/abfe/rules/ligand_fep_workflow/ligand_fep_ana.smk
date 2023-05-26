from abfe.scripts.free_energy import analysis
from abfe.utils import tools, mdp
import os
import shutil

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
        boresch_dat = run_path+"/complex/equil-mdsim/boreschcalc/dG_off.dat",
    params:
        ana_loc=run_path+"/ligand/fep/ana",
    output:
        ligand_json=run_path+"/ligand/fep/ana/dg_ligand_contributions.json"
    threads: threads
    run:
        # Make directory
        tools.makedirs(params.ana_loc)        
        # Get the simulaiton temperature from the prod.mdp of the state 0 of vdw
        temperature = float(mdp.MDP().from_file(input.mdp_vdw_0_prod).parameters['ref-t'].split()[0]) 
        analysis.get_dG_contributions(
            boresch_data = input.boresch_dat,
            out_json_path = output.ligand_json,
            # Check if it is necessary to remove some initial burning simulation time
            lower = None,
            upper = None,
            min_samples = 500,
            temperature = temperature,
            vdw = sorted(input.xvg_vdw_loc, key=lambda x: int(os.path.normpath(x).split(os.path.sep)[-3].split('.')[-1])),
            coul = sorted(input.xvg_coul_loc, key=lambda x: int(os.path.normpath(x).split(os.path.sep)[-3].split('.')[-1])),
        )
