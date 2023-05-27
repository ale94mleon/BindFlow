from abfe.scripts.free_energy import analysis
run_path = config["run_path"]


rule fep_get_dg_cycle:
    input:
        complex_json=run_path+"/complex/fep/ana/dg_complex_contributions.json",
        ligand_json=run_path+"/ligand/fep/ana/dg_ligand_contributions.json",
    output:
        out_file_path=run_path+"/dG_results.csv",
    run:
        analysis.get_dg_cycle(
            ligand_contributions = input.ligand_json,
            complex_contributions = input.complex_json,
            out_csv = output.out_file_path
        )

