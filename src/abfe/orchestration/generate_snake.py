from abfe import rules


def get_eq_res():
    equil_snake = rules.ligand_equilibration_workflow_path

    cmd = [
        "#Do Equilibration",
        "include: \'" + equil_snake + "\'",
        "rule check_equilib:",
        "   input:",
        "       gro_complex=run_path+\"/complex/equil-mdsim/boreschcalc/ClosestRestraintFrame.gro\",",
        "       top_complex=run_path+\"/complex/equil-mdsim/boreschcalc/BoreschRestraint.top\",",
        "       dG_off_complex=run_path+\"/complex/equil-mdsim/boreschcalc/dG_off.dat\",",
        "       gro_ligand=run_path+\"/ligand/equil-mdsim/npt_equil2/npt_equil2.gro\",",
        "       cpt_ligand=run_path+\"/ligand/equil-mdsim/npt_equil2/npt_equil2.cpt\"",
        ""
    ]
    return "\n".join(cmd)


def get_fep_res() -> str:
    fep_snake = rules.ligand_fep_workflow_path

    cmd = [
        "#Do FEP",
        "include: \'" + fep_snake + "\'",
        "rule all:",
        "    input:",
        "        dG_path=run_path+\"/dG_results.csv\"",
        ""
    ]
    return "\n".join(cmd)


def get_all_eq_fep_res() -> str:
    cmd = [
        "#DO",
        "## Do Equilibration",
        "include: \'" + rules.ligand_equilibration_workflow_path + "\'",
        "",
        "## Do FEP",
        "include: \'" + rules.ligand_fep_workflow_path + "\'",
        "",
        "## Do Check all results",
        "rule abfe_ligand_result:",
        "    input:",
        "        dG_path=run_path+\"/dG_results.csv\"",
        ""
    ]
    return "\n".join(cmd)


def get_superFlow() -> str:
    cmd = [
        "#DO",
        "## Start the launching process",
        "include: \'" + rules.receptor_workflow_path + "\'",
        "",
        "## Do Check all results",
        "rule abfe_recptor_result:",
        "    input:",
        "        dG_path=approach_path+\"/abfe_results.csv\"",
        ""
    ]
    return "\n".join(cmd)


def generate_approach_snake_file(out_file_path: str,
                                 conf_file_path: str):
    load_conf_file = "configfile: " \
                     "\'" + conf_file_path + "\'\napproach_path =config['out_approach_path']\n"

    full_job = get_superFlow()

    file_str = "\n".join(["#Load Config:", load_conf_file, full_job])

    out_file_IO = open(out_file_path, "w")
    out_file_IO.write(file_str)
    out_file_IO.close()


def generate_snake_file(out_file_path: str,
                        conf_file_path: str):
    load_conf_file = "configfile: \'" + conf_file_path + "\'\nrun_path = config['run_path']\n"

    full_job = get_all_eq_fep_res()

    file_str = "\n".join(["#Load Config:", load_conf_file,
                          full_job,
                          ])

    out_file_IO = open(out_file_path, "w")
    out_file_IO.write(file_str)
    out_file_IO.close()
