from bindflow.free_energy import gather_results

# Gather Results
rule gather_receptor_results:
    input:
        mmxbsa_csv=expand(out_approach_path+"/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/mmxbsa.csv", ligand_name = config['ligand_names'], replica = list(map(str, range(1,1 + config['replicas']))), sample = list(map(str, range(1,1 + config["samples"]))))
    output:
        out_dg_file=out_approach_path+"/mmxbsa_results.csv",
        out_raw_file=out_approach_path+"/mmxbsa_results_raw.csv"
    run:
        full_df = gather_results.get_raw_mmxbsa_dgs(
            root_folder_path=out_approach_path,
            out_csv=output.out_raw_file
        )
        gather_results.get_all_mmxbsa_dgs(
            full_df=full_df,
            columns_to_process=None,
            out_csv=output.out_dg_file
        )