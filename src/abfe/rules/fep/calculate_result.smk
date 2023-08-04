from abfe.free_energy import analysis
approach_path = config["out_approach_path"]

rule fep_get_dg_cycle:
    input:
        complex_json=approach_path + "/{ligand_name}/{replica}/complex/fep/ana/dg_complex_contributions.json",
        ligand_json=approach_path + "/{ligand_name}/{replica}/ligand/fep/ana/dg_ligand_contributions.json",
    output:
        out_file_path=approach_path + "/{ligand_name}/{replica}/dG_results.csv",
    run:
        analysis.get_dg_cycle(
            ligand_contributions = input.ligand_json,
            complex_contributions = input.complex_json,
            out_csv = output.out_file_path
        )