from bindflow.free_energy import analysis

approach_path = config["out_approach_path"]

rule mmpbsa_get_dg_cycle:
    input:
        complex_equi_finished = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.finished"
    output:
        out_file_path = approach_path + "/{ligand_name}/{replica}/dG_results.csv",
    params:
        in_gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.gro",
        in_cpt = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.cpt",
        in_tpr = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.tpr",
        in_xtc = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.xtc",
    run:
        # TODO Implement the function to get the free energy of binding
        pass