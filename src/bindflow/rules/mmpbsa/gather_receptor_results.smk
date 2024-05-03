from bindflow.free_energy import gather_results

out_approach_path = config["out_approach_path"]

# Gather Results
rule gather_receptor_results:
    input:
        mmxbsa_csv = expand(approach_path + "/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/mmxbsa.csv", ligand_name = config['ligand_names'], replica = list(map(str, range(1,1 + config['replicas']))), sample = list(map(str, range(1,1 + config["samples"]))))
    output:
        out_dg_file=out_approach_path + "/mmxbsa_results.csv",
        out_raw_file=out_approach_path + "/mmxbsa_results_raw.csv",
    run:
        pass
        # Implementation