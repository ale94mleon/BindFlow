#Final Check Job
approach_path = config["out_approach_path"]

rule calculate_ABFE:
    input:
        in_ligand_gro= approach_path + "/{ligand_name}/input/ligand/ligand.gro",
        in_ligand_top= approach_path + "/{ligand_name}/input/ligand/ligand.top",
        in_complex_gro = approach_path + "/{ligand_name}/input/complex/complex.gro",
        in_complex_top = approach_path + "/{ligand_name}/input/complex/complex.top",
        output_dir = approach_path + "/{ligand_name}/{replica}",

    output:
        out_dg = approach_path + "/{ligand_name}/{replica}/dG_results.csv"

    shell:
        """
            cd {input.output_dir}
            ./job.sh
        """