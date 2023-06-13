from abfe.utils.tools import gmx_runner

# Common to all the sub-workflows ligand/replica
input_path = config['input_data_path']
run_path = config["run_path"]
threads = config['threads']
retries = config['retries']
load_dependencies = config['extra_directives']['dependencies']
mdrun_extra = config['extra_directives']['mdrun']

rule equil_ligand_00_min:
    input:
        top=input_path+"/ligand/ligand.top",
        gro=input_path+"/ligand/ligand.gro",
        mdp=run_path+"/ligand/equil-mdsim/00_min/00_min.mdp"
    params:
        run_dir=run_path+"/ligand/equil-mdsim/00_min"
    output:
        gro=run_path+"/ligand/equil-mdsim/00_min/00_min.gro"
    threads: threads
    retries: retries
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

rule equil_ligand_01_nvt:
    input:
        top=input_path+"/ligand/ligand.top",
        gro=run_path+"/ligand/equil-mdsim/00_min/00_min.gro",
        mdp=run_path+"/ligand/equil-mdsim/01_nvt/01_nvt.mdp",
    params:
        run_dir=run_path+"/ligand/equil-mdsim/01_nvt"
    output:
        gro=run_path+"/ligand/equil-mdsim/01_nvt/01_nvt.gro",
        cpt=run_path+"/ligand/equil-mdsim/01_nvt/01_nvt.cpt"
    threads: threads
    retries: retries
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

rule equil_ligand_02_nvt:
    input:
        top=input_path+"/ligand/ligand.top",
        gro=run_path+"/ligand/equil-mdsim/01_nvt/01_nvt.gro",
        cpt=run_path+"/ligand/equil-mdsim/01_nvt/01_nvt.cpt",
        mdp=run_path+"/ligand/equil-mdsim/02_nvt/02_nvt.mdp",
    params:
        run_dir=run_path+"/ligand/equil-mdsim/02_nvt"
    output:
        gro=run_path+"/ligand/equil-mdsim/02_nvt/02_nvt.gro",
        cpt=run_path+"/ligand/equil-mdsim/02_nvt/02_nvt.cpt"
    threads: threads
    retries: retries
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

rule equil_ligand_03_npt:
    input:
        top=input_path+"/ligand/ligand.top",
        gro=run_path+"/ligand/equil-mdsim/02_nvt/02_nvt.gro",
        cpt=run_path+"/ligand/equil-mdsim/02_nvt/02_nvt.cpt",
        mdp=run_path+"/ligand/equil-mdsim/03_npt/03_npt.mdp",
    params:
        run_dir=run_path+"/ligand/equil-mdsim/03_npt"
    output:
        gro=run_path+"/ligand/equil-mdsim/03_npt/03_npt.gro",
        cpt=run_path+"/ligand/equil-mdsim/03_npt/03_npt.cpt"
    threads: threads
    retries: retries
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

rule equil_ligand_04_npt:
    input:
        top=input_path+"/ligand/ligand.top",
        gro=run_path+"/ligand/equil-mdsim/03_npt/03_npt.gro",
        cpt=run_path+"/ligand/equil-mdsim/03_npt/03_npt.cpt",
        mdp=run_path+"/ligand/equil-mdsim/04_npt/04_npt.mdp",
    params:
        run_dir=run_path+"/ligand/equil-mdsim/04_npt"
    output:
        gro=run_path+"/ligand/equil-mdsim/04_npt/04_npt.gro",
        cpt=run_path+"/ligand/equil-mdsim/04_npt/04_npt.cpt",
    threads: threads
    retries: retries
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
