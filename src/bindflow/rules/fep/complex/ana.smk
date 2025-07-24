from bindflow.free_energy import fep_analysis
# Ana
rule fep_ana_get_dg_complex_contributions:
    input:
        # Make sure that the simulation ends properly
        finished_vdw=expand(out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/vdw.{state}/prod/prod.finished", state=range(len(config['lambdas']['complex']['vdw'])), allow_missing=True),
        finished_coul=expand(out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/coul.{state}/prod/prod.finished", state=range(len(config['lambdas']['complex']['coul'])), allow_missing=True),
        finished_bonded=expand(out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/bonded.{state}/prod/prod.finished", state=range(len(config['lambdas']['complex']['bonded'])), allow_missing=True),
        # Boresch correction
        boresch_dat=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/boreschcalc/dG_off.dat",
        # To get the simulation temperature
        mdp=expand(out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{sim_type}.{state}/prod/prod.mdp", state=range(len(config['lambdas']['complex']['bonded'])), sim_type=['vdw', 'coul', 'bonded'], allow_missing=True)
    params:
        xvg_vdw=expand(out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/vdw.{state}/prod/prod.xvg", state=range(len(config['lambdas']['complex']['vdw'])), allow_missing=True),
        xvg_coul=expand(out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/coul.{state}/prod/prod.xvg", state=range(len(config['lambdas']['complex']['coul'])), allow_missing=True),
        xvg_bonded=expand(out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/bonded.{state}/prod/prod.xvg", state=range(len(config['lambdas']['complex']['bonded'])), allow_missing=True),
        ana=out_approach_path+"/{ligand_name}/{replica}/complex/fep/ana",
    output:
        complex_json=out_approach_path+"/{ligand_name}/{replica}/complex/fep/ana/dg_complex_contributions.json"
    threads: threads # TODO: Sometimes the rule hang for a long time
    run:
        # Make directory
        Path(params.ana).mkdir(exist_ok=True, parents=True)
        # Get simulation temperature from any prod.mdp file (all should have the same)
        mdp_params = mdp.MDP().from_file(input.mdp[0]).parameters
        if 'ref-t' in mdp_params:
            temperature = float(mdp_params['ref-t'].split()[0])
        elif 'ref_t' in mdp_params:
            temperature = float(mdp_params['ref_t'].split()[0])
        fep_analysis.get_dG_contributions(
            boresch_data=input.boresch_dat,
            out_json_path=output.complex_json,
            # Check if it is necessary to remove some initial burning simulation time
            lower=None,
            upper=None,
            min_samples=500,
            temperature=temperature,
            # convergency_plots_prefix = params.ana + "/complex_",
            convergency_plots_prefix=None,
            # Sort the paths
            vdw=sorted(params.xvg_vdw, key=lambda x: int(Path(x).parts[-3].split('.')[-1])),
            coul=sorted(params.xvg_coul, key=lambda x: int(Path(x).parts[-3].split('.')[-1])),
            bonded=sorted(params.xvg_bonded, key=lambda x: int(Path(x).parts[-3].split('.')[-1])),
        )

