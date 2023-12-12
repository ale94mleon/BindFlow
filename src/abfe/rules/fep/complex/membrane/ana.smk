from abfe.free_energy import analysis
from abfe.utils import tools
from abfe.mdp import mdp
import os

threads = config['threads']

# Ana
rule fep_ana_get_dg_complex_contributions:
    input:
        # Make sure that the simualtion ends properly
        finished_vdw_loc = expand(approach_path + "/{ligand_name}/{replica}/complex/fep/simulation/vdw.{state}/prod/prod.finished", state = range(len(config['lambdas']['complex']['vdw'])), allow_missing = True),
        finished_coul_loc = expand(approach_path + "/{ligand_name}/{replica}/complex/fep/simulation/coul.{state}/prod/prod.finished", state = range(len(config['lambdas']['complex']['coul'])), allow_missing = True),
        finished_bonded_loc = expand(approach_path + "/{ligand_name}/{replica}/complex/fep/simulation/bonded.{state}/prod/prod.finished", state = range(len(config['lambdas']['complex']['bonded'])), allow_missing = True),
        # Boresch correction
        boresch_dat = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/boreschcalc/dG_off.dat",
        # To get the simulaiton temperature
        mdp_vdw_0_prod = approach_path + "/{ligand_name}/{replica}/complex/fep/simulation/vdw.0/prod/prod.mdp",
    params:
        xvg_vdw_loc = expand(approach_path + "/{ligand_name}/{replica}/complex/fep/simulation/vdw.{state}/prod/prod.xvg", state = range(len(config['lambdas']['complex']['vdw'])), allow_missing = True),
        xvg_coul_loc = expand(approach_path + "/{ligand_name}/{replica}/complex/fep/simulation/coul.{state}/prod/prod.xvg", state = range(len(config['lambdas']['complex']['coul'])), allow_missing = True),
        xvg_bonded_loc = expand(approach_path + "/{ligand_name}/{replica}/complex/fep/simulation/bonded.{state}/prod/prod.xvg", state = range(len(config['lambdas']['complex']['bonded'])), allow_missing = True),
        ana_loc = approach_path + "/{ligand_name}/{replica}/complex/fep/ana",
    output:
        complex_json = approach_path + "/{ligand_name}/{replica}/complex/fep/ana/dg_complex_contributions.json"
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
            out_json_path = output.complex_json,
            # Check if it is necessary to remove some initial burning simulation time
            lower = None,
            upper = None,
            min_samples = 500,
            temperature = temperature,
            # convergency_plots_prefix = params.ana_loc + "/complex_",
            convergency_plots_prefix = None,
            # Sort the paths
            vdw = sorted(params.xvg_vdw_loc, key = lambda x: int(os.path.normpath(x).split(os.path.sep)[-3].split('.')[-1])),
            coul = sorted(params.xvg_coul_loc, key = lambda x: int(os.path.normpath(x).split(os.path.sep)[-3].split('.')[-1])),
            bonded = sorted(params.xvg_bonded_loc, key = lambda x: int(os.path.normpath(x).split(os.path.sep)[-3].split('.')[-1])),
        )

