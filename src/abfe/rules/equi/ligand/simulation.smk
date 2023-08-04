from abfe.utils import tools

# Common to all the sub-workflows ligand/replica
approach_path = config["out_approach_path"]
threads = config['threads']
retries = config['retries']
load_dependencies = config['extra_directives']['dependencies']
mdrun_extra = config['extra_directives']['mdrun']

rule equil_ligand_00_min:
    input:
        top=approach_path + "/{ligand_name}/input/ligand/ligand.top",
        gro=approach_path + "/{ligand_name}/input/ligand/ligand.gro",
        mdp=approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/00_min/00_min.mdp"
    params:
        run_dir=approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/00_min"
    output:
        gro=approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/00_min/00_min.gro"
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
        top=approach_path + "/{ligand_name}/input/ligand/ligand.top",
        gro=approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/00_min/00_min.gro",
        mdp=approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/01_nvt/01_nvt.mdp",
    params:
        output={
            'gro':approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/01_nvt/01_nvt.gro",
            'cpt':approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/01_nvt/01_nvt.cpt"
        },
        run_dir=approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/01_nvt"
    output:
        finished=approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/01_nvt/01_nvt.finished",
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
        top=approach_path + "/{ligand_name}/input/ligand/ligand.top",
        finished=approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/01_nvt/01_nvt.finished",
        mdp=approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/02_nvt/02_nvt.mdp",
    params:
        input={
            'gro':approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/01_nvt/01_nvt.gro",
            'cpt':approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/01_nvt/01_nvt.cpt"
        },
        output={
            'gro':approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/02_nvt/02_nvt.gro",
            'cpt':approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/02_nvt/02_nvt.cpt"
        },
        run_dir=approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/02_nvt"
    output:
        finished=approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/02_nvt/02_nvt.finished",
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
        top=approach_path + "/{ligand_name}/input/ligand/ligand.top",
        finished=approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/02_nvt/02_nvt.finished",
        mdp=approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/03_npt/03_npt.mdp",
    params:
        input={
            'gro':approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/02_nvt/02_nvt.gro",
            'cpt':approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/02_nvt/02_nvt.cpt"
        },
        output={
            'gro':approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/03_npt/03_npt.gro",
            'cpt':approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/03_npt/03_npt.cpt"
        },
        run_dir=approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/03_npt"
    output:
        finished=approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/03_npt/03_npt.finished",
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
        top=approach_path + "/{ligand_name}/input/ligand/ligand.top",
        finished=approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/03_npt/03_npt.finished",
        mdp=approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/prod/prod.mdp",
    params:
        input={
            'gro':approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/03_npt/03_npt.gro",
            'cpt':approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/03_npt/03_npt.cpt"
        },
        output={
            'gro':approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/prod/prod.gro",
            'cpt':approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/prod/prod.cpt",
        },
        run_dir=approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/prod"
    output:
        finished=approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/prod/prod.finished",
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