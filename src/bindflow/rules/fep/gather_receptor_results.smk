from bindflow.free_energy import gather_results

out_approach_path = config["out_approach_path"]

# Gather Results
rule gather_receptor_results:
    input:
        prior_result_paths = expand(out_approach_path + "/{ligand_names}/{replica}/dG_results.csv",ligand_names = config['ligand_names'],replica = list(map(str, range(1,1 + config['replicas']))))
    output:
        out_dg_file=out_approach_path + "/abfe_results.csv",
        out_raw_file=out_approach_path + "/abfe_results_raw.csv",
    run:
        gather_results.get_all_dgs(
            root_folder_path=out_approach_path,
            out_csv=os.path.join(output.out_dg_file)
        )
        gather_results.get_raw_data(
            root_folder_path=out_approach_path,
            out_csv=os.path.join(output.out_raw_file)
        )