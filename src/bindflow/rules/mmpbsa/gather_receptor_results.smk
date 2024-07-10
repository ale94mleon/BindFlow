from bindflow.free_energy import mmxbsa_analysis
import pandas as pd

# Gather Results
rule gather_receptor_results:
    input:
        mmxbsa_csv=expand(out_approach_path+"/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/mmxbsa.csv", ligand_name = config['ligand_names'], replica = list(map(str, range(1,1 + config['replicas']))), sample = list(map(str, range(1,1 + config["samples"]))))
    output:
        out_dg_file=out_approach_path+"/mmxbsa_results.csv",
        out_raw_file=out_approach_path+"/mmxbsa_results_raw.csv"
    run:
        collected_dfs = []
        for inp_file in input.mmxbsa_csv: # collecting all of the mmxbsa.csv files and concatenating into 1 file
            string_data = inp_file.removeprefix(out_approach_path).split("/")
            ligand_name, replica, sample = string_data[1], string_data[2], string_data[6].removeprefix("rep.")
            collected_dfs.append(mmxbsa_analysis.convert_format_flatten(pd.read_csv(inp_file), ligand_name, replica, sample))
        full_df = pd.concat(collected_dfs, ignore_index=True)
        full_df.to_csv(output.out_raw_file, index=False)
        
        mmxbsa_analysis.prettify_df(full_df).to_csv(output.out_dg_file, index=False)