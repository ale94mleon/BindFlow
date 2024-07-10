from bindflow.free_energy import analysis
# Ana

# TODO Here is the main issue!!!! How can I isolate for each ligand and each replica
# If I use them as widlcards all the simulaiton must end, and that is not what we want
rule fep_ana_get_dg_ligand_contributions:
    input:
        # Make sure that the simualtion ends properly
        finished_vdw_loc=expand(out_approach_path+"/{ligand_name}/{replica}/ligand/fep/simulation/vdw.{state}/prod/prod.finished", state=range(len(config['lambdas']['ligand']['vdw'])), allow_missing=True),
        finished_coul_loc=expand(out_approach_path+"/{ligand_name}/{replica}/ligand/fep/simulation/coul.{state}/prod/prod.finished", state=range(len(config['lambdas']['ligand']['coul'])), allow_missing=True),
        # To get the simulaiton temperature
        mdp_vdw_0_prod=out_approach_path+"/{ligand_name}/{replica}/ligand/fep/simulation/vdw.0/prod/prod.mdp",
    params:
        #  TODO finished_vdw_loc is needed to connect the rule dependencies, but xvg_vdw_loc is the thing that I need and they could also be passed as input. if finished is there xvg should also be there.
        xvg_vdw_loc=expand(out_approach_path+"/{ligand_name}/{replica}/ligand/fep/simulation/vdw.{state}/prod/prod.xvg", state=range(len(config['lambdas']['ligand']['vdw'])), allow_missing=True),
        xvg_coul_loc=expand(out_approach_path+"/{ligand_name}/{replica}/ligand/fep/simulation/coul.{state}/prod/prod.xvg", state=range(len(config['lambdas']['ligand']['coul'])), allow_missing=True),
        ana_loc=out_approach_path+"/{ligand_name}/{replica}/ligand/fep/ana",
    output:
        ligand_json=out_approach_path+"/{ligand_name}/{replica}/ligand/fep/ana/dg_ligand_contributions.json"
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
            boresch_data=None,
            out_json_path=output.ligand_json,
            # Check if it is necessary to remove some initial burning simulation time
            lower=None,
            upper=None,
            min_samples=500,
            temperature=temperature,
            # convergency_plots_prefix = params.ana_loc + "/ligand_",
            convergency_plots_prefix=None,
            # Sort the paths
            vdw=sorted(params.xvg_vdw_loc, key=lambda x: int(os.path.normpath(x).split(os.path.sep)[-3].split('.')[-1])),
            coul=sorted(params.xvg_coul_loc, key=lambda x: int(os.path.normpath(x).split(os.path.sep)[-3].split('.')[-1])),
        )