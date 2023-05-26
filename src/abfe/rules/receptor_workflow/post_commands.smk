from abfe import scripts

#Final Check Job
rule gather_receptor_results:
    input:
        approach_path=config["out_approach_path"],
        prior_result_paths=expand(config["out_approach_path"] + "/{ligand_names}/{replica}/dG_results.csv",ligand_names=config['ligand_names'],replica=list(map(str, range(1,1 + config['replicas']))), allow_missing=True)
    params:
        script_dir = scripts.root_path
    output:
        out_dG_File=config["out_approach_path"]+"/abfe_results.csv",
    # TODO, I have to modify this part now
    shell:
        "python {params.script_dir}/final_receptor_results.py --in_root_dir {input.approach_path} --out_dir {input.approach_path}"