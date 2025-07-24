from bindflow.free_energy import fep_analysis
# Ana

# TODO Here is the main issue!!!! How can I isolate for each ligand and each replica
# If I use them as widlcards all the simulaiton must end, and that is not what we want
rule fep_ana_get_dg_ligand_contributions:
    input:
        # Make sure that the simualtion ends properly
        finished_vdw=expand(out_approach_path+"/{ligand_name}/{replica}/ligand/fep/simulation/vdw.{state}/prod/prod.finished", state=range(len(config['lambdas']['ligand']['vdw'])), allow_missing=True),
        finished_coul=expand(out_approach_path+"/{ligand_name}/{replica}/ligand/fep/simulation/coul.{state}/prod/prod.finished", state=range(len(config['lambdas']['ligand']['coul'])), allow_missing=True),
        # To get the simulation temperature
        mdp=expand(out_approach_path+"/{ligand_name}/{replica}/ligand/fep/simulation/{sim_type}.{state}/prod/prod.mdp", state=range(len(config['lambdas']['complex']['bonded'])), sim_type=['vdw', 'coul'], allow_missing=True)
    params:
        #  TODO finished_vdw is needed to connect the rule dependencies, but xvg_vdw is the thing that I need and they could also be passed as input. if finished is there xvg should also be there.
        xvg_vdw=expand(out_approach_path+"/{ligand_name}/{replica}/ligand/fep/simulation/vdw.{state}/prod/prod.xvg", state=range(len(config['lambdas']['ligand']['vdw'])), allow_missing=True),
        xvg_coul=expand(out_approach_path+"/{ligand_name}/{replica}/ligand/fep/simulation/coul.{state}/prod/prod.xvg", state=range(len(config['lambdas']['ligand']['coul'])), allow_missing=True),
        ana=out_approach_path+"/{ligand_name}/{replica}/ligand/fep/ana",
    output:
        ligand_json=out_approach_path+"/{ligand_name}/{replica}/ligand/fep/ana/dg_ligand_contributions.json"
    threads: threads # TODO: Sometimes the rule hang for a long time
    run:
        # Make directory
        # Get simulation temperature from any prod.mdp file (all should have the same)
        mdp_params = mdp.MDP().from_file(input.mdp[0]).parameters
        if 'ref-t' in mdp_params:
            temperature = float(mdp_params['ref-t'].split()[0])
        elif 'ref_t' in mdp_params:
            temperature = float(mdp_params['ref_t'].split()[0])

        fep_analysis.get_dG_contributions(
            boresch_data=None,
            out_json_path=output.ligand_json,
            # Check if it is necessary to remove some initial burning simulation time
            lower=None,
            upper=None,
            min_samples=500,
            temperature=temperature,
            # convergency_plots_prefix = params.ana + "/ligand_",
            convergency_plots_prefix=None,
            # Sort the paths
            vdw=sorted(params.xvg_vdw, key=lambda x: int(Path(x).parts[-3].split('.')[-1])),
            coul=sorted(params.xvg_coul, key=lambda x: int(Path(x).parts[-3].split('.')[-1])),
        )