from abfe import scripts

#Final Check Job
approach_path = config["out_approach_path"]
ligand_names = config['ligand_names']
replicas =  config['replicas']
replica_list = list(map(str, range(1,1 + replicas)))

rule gather_receptor_results:
    input:
        approach_path=approach_path,
        prior_result_paths=expand(approach_path + "/{ligand_names}/{replica}/dG_results.csv",ligand_names=ligand_names,replica=replica_list, allow_missing=True)
    params:
        script_dir = scripts.root_path
    output:
        out_dG_File=approach_path+"/abfe_results.csv",
    shell:
        "python {params.script_dir}/final_receptor_results.py --in_root_dir {input.approach_path} --out_dir {input.approach_path}"