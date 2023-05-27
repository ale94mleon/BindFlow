from abfe.utils.tools import gmx_runner

# Common to all the sub-workflows ligand/replica
input_path = config['input_data_path']
run_path = config["run_path"]
threads = config['threads']
num_retries = config['num_retries']
load_dependencies = config['extra_directives']['dependencies']
mdrun_extra = config['extra_directives']['mdrun']

rule equil_ligand_emin:
    input:
        top=input_path+"/ligand/ligand.top",
        gro=input_path+"/ligand/ligand.gro",
        mdp=run_path+"/ligand/equil-mdsim/emin/emin.mdp"
    params:
        run_dir=run_path+"/ligand/equil-mdsim/emin"
    output:
        gro=run_path+"/ligand/equil-mdsim/emin/emin.gro"
    threads: threads
    retries: num_retries
    run:
        gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = input.gro,
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['ligand']
        )

rule equil_ligand_nvt_heat:
    input:
        top=input_path+"/ligand/ligand.top",
        gro=run_path+"/ligand/equil-mdsim/emin/emin.gro",
        mdp=run_path+"/ligand/equil-mdsim/nvt_heat/nvt_heat.mdp",
    params:
        run_dir=run_path+"/ligand/equil-mdsim/nvt_heat"
    output:
        gro=run_path+"/ligand/equil-mdsim/nvt_heat/nvt_heat.gro",
        cpt=run_path+"/ligand/equil-mdsim/nvt_heat/nvt_heat.cpt"
    threads: threads
    retries: num_retries
    run:
        gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = input.gro,
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['ligand']
        )

rule equil_ligand_npt_eq1:
    input:
        top=input_path+"/ligand/ligand.top",
        gro=run_path+"/ligand/equil-mdsim/nvt_heat/nvt_heat.gro",
        cpt=run_path+"/ligand/equil-mdsim/nvt_heat/nvt_heat.cpt",
        mdp=run_path+"/ligand/equil-mdsim/npt_equil1/npt_equil1.mdp",
    params:
        run_dir=run_path+"/ligand/equil-mdsim/npt_equil1"
    output:
        gro=run_path+"/ligand/equil-mdsim/npt_equil1/npt_equil1.gro",
        cpt=run_path+"/ligand/equil-mdsim/npt_equil1/npt_equil1.cpt"
    threads: threads
    retries: num_retries
    run:
        gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = input.gro,
            checkpoint = input.cpt,
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['ligand']
        )

rule equil_ligand_npt_eq2:
    input:
        top=input_path+"/ligand/ligand.top",
        gro=run_path+"/ligand/equil-mdsim/npt_equil1/npt_equil1.gro",
        cpt=run_path+"/ligand/equil-mdsim/npt_equil1/npt_equil1.cpt",
        mdp=run_path+"/ligand/equil-mdsim/npt_equil2/npt_equil2.mdp",
    params:
        run_dir=run_path+"/ligand/equil-mdsim/npt_equil2"
    output:
        gro=run_path+"/ligand/equil-mdsim/npt_equil2/npt_equil2.gro",
        cpt=run_path+"/ligand/equil-mdsim/npt_equil2/npt_equil2.cpt",
    threads: threads
    retries: num_retries
    run:
        gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = input.gro,
            checkpoint = input.cpt,
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['ligand']
        )
