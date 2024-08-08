from bindflow.free_energy import fep_analysis

rule fep_get_dg_cycle:
    input:
        complex_json=out_approach_path+"/{ligand_name}/{replica}/complex/fep/ana/dg_complex_contributions.json",
        ligand_json=out_approach_path+"/{ligand_name}/{replica}/ligand/fep/ana/dg_ligand_contributions.json",
    output:
        out_file_path=out_approach_path+"/{ligand_name}/{replica}/dG_results.csv",
    run:
        fep_analysis.get_dg_cycle(
            ligand_contributions=input.ligand_json,
            complex_contributions=input.complex_json,
            out_csv=output.out_file_path
        )