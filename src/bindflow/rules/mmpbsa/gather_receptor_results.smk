from bindflow.free_energy import gather_results

out_approach_path = config["out_approach_path"]

# Gather Results
rule gather_receptor_results:
    input:
        prior_result_paths = expand(out_approach_path + "/{ligand_names}/{replica}/dG_mmxbsa_results.csv",ligand_names = config['ligand_names'],replica = list(map(str, range(1,1 + config['replicas']))))
    output:
        out_dg_file=out_approach_path + "/mmxbsa_results.csv",
        out_raw_file=out_approach_path + "/mmxbsa_results_raw.csv",
    run:
        pass
        # Implementation