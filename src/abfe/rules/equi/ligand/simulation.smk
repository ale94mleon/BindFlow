from abfe.utils import tools

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
        tools.gmx_runner(
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
        output={
            'gro':run_path+"/ligand/equil-mdsim/01_nvt/01_nvt.gro",
            'cpt':run_path+"/ligand/equil-mdsim/01_nvt/01_nvt.cpt"
        },
        run_dir=run_path+"/ligand/equil-mdsim/01_nvt"
    output:
        finished=run_path+"/ligand/equil-mdsim/01_nvt/01_nvt.finished",
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = input.gro,
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['ligand']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths = params.output.values(), raise_error = True, out = output.finished)

rule equil_ligand_02_nvt:
    input:
        top=input_path+"/ligand/ligand.top",
        finished=run_path+"/ligand/equil-mdsim/01_nvt/01_nvt.finished",
        mdp=run_path+"/ligand/equil-mdsim/02_nvt/02_nvt.mdp",
    params:
        input={
            'gro':run_path+"/ligand/equil-mdsim/01_nvt/01_nvt.gro",
            'cpt':run_path+"/ligand/equil-mdsim/01_nvt/01_nvt.cpt"
        },
        output={
            'gro':run_path+"/ligand/equil-mdsim/02_nvt/02_nvt.gro",
            'cpt':run_path+"/ligand/equil-mdsim/02_nvt/02_nvt.cpt"
        },
        run_dir=run_path+"/ligand/equil-mdsim/02_nvt"
    output:
        finished=run_path+"/ligand/equil-mdsim/02_nvt/02_nvt.finished",
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = params.input['gro'],
            checkpoint = params.input['cpt'],
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['ligand']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths = params.output.values(), raise_error = True, out = output.finished)

rule equil_ligand_03_npt:
    input:
        top=input_path+"/ligand/ligand.top",
        finished=run_path+"/ligand/equil-mdsim/02_nvt/02_nvt.finished",
        mdp=run_path+"/ligand/equil-mdsim/03_npt/03_npt.mdp",
    params:
        input={
            'gro':run_path+"/ligand/equil-mdsim/02_nvt/02_nvt.gro",
            'cpt':run_path+"/ligand/equil-mdsim/02_nvt/02_nvt.cpt"
        },
        output={
            'gro':run_path+"/ligand/equil-mdsim/03_npt/03_npt.gro",
            'cpt':run_path+"/ligand/equil-mdsim/03_npt/03_npt.cpt"
        },
        run_dir=run_path+"/ligand/equil-mdsim/03_npt"
    output:
        finished=run_path+"/ligand/equil-mdsim/03_npt/03_npt.finished",
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = params.input['gro'],
            checkpoint = params.input['cpt'],
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['ligand']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths = params.output.values(), raise_error = True, out = output.finished)

rule equil_ligand_prod:
    input:
        top=input_path+"/ligand/ligand.top",
        finished=run_path+"/ligand/equil-mdsim/03_npt/03_npt.finished",
        mdp=run_path+"/ligand/equil-mdsim/prod/prod.mdp",
    params:
        input={
            'gro':run_path+"/ligand/equil-mdsim/03_npt/03_npt.gro",
            'cpt':run_path+"/ligand/equil-mdsim/03_npt/03_npt.cpt"
        },
        output={
            'gro':run_path+"/ligand/equil-mdsim/prod/prod.gro",
            'cpt':run_path+"/ligand/equil-mdsim/prod/prod.cpt",
        },
        run_dir=run_path+"/ligand/equil-mdsim/prod"
    output:
        finished=run_path+"/ligand/equil-mdsim/prod/prod.finished",
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = params.input['gro'],
            checkpoint = params.input['cpt'],
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['ligand']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths = params.output.values(), raise_error = True, out = output.finished)