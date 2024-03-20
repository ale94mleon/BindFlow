from bindflow.utils import tools

# Common to all the sub-workflows ligand/replica
approach_path = config["out_approach_path"]
threads = config['threads']
retries = config['retries']
load_dependencies = config['extra_directives']['dependencies']
mdrun_extra = config['extra_directives']['mdrun']

rule equil_complex_00_min:
    input:
        top = approach_path + "/{ligand_name}/input/complex/complex.top",
        ndx = approach_path + "/{ligand_name}/input/complex/index.ndx",
        gro = approach_path + "/{ligand_name}/input/complex/complex.gro",
        mdp = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/00_min/00_min.mdp"
    params:
        run_dir = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/00_min"
    output:
        gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/00_min/00_min.gro"
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = input.gro,
            index = input.ndx,
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['complex']
        )

rule equil_complex_01_nvt:
    input:
        top = approach_path + "/{ligand_name}/input/complex/complex.top",
        ndx = approach_path + "/{ligand_name}/input/complex/index.ndx",
        gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/00_min/00_min.gro",
        mdp = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt/01_nvt.mdp",
    params:
        out_gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt/01_nvt.gro",
        out_cpt = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt/01_nvt.cpt",
        run_dir = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt",
    output:
        finished = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt/01_nvt.finished",
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = input.gro,
            index = input.ndx,
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['complex']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths = [params.out_gro, params.out_cpt], raise_error = True, out = output.finished)

rule equil_complex_02_npt:
    input:
        top = approach_path + "/{ligand_name}/input/complex/complex.top",
        ndx = approach_path + "/{ligand_name}/input/complex/index.ndx",
        finished = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt/01_nvt.finished",
        mdp = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/02_npt/02_npt.mdp",
    params:
        in_gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt/01_nvt.gro",
        in_cpt = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt/01_nvt.cpt",
        out_gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/02_npt/02_npt.gro",
        out_cpt = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/02_npt/02_npt.cpt",
        run_dir = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/02_npt",
    output:
        finished = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/02_npt/02_npt.finished",
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = params.in_gro,
            index = input.ndx,
            checkpoint = params.in_cpt,
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['complex']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths = [params.out_gro, params.out_cpt], raise_error = True, out = output.finished)

rule equil_complex_03_npt:
    input:
        top = approach_path + "/{ligand_name}/input/complex/complex.top",
        ndx = approach_path + "/{ligand_name}/input/complex/index.ndx",
        finished = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/02_npt/02_npt.finished",
        mdp = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/03_npt/03_npt.mdp",
    params:
        in_gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/02_npt/02_npt.gro",
        in_cpt = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/02_npt/02_npt.cpt",
        out_gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/03_npt/03_npt.gro",
        out_cpt = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/03_npt/03_npt.cpt",
        run_dir = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/03_npt",
    output:
        finished = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/03_npt/03_npt.finished",
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = params.in_gro,
            index = input.ndx,
            checkpoint = params.in_cpt,
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['complex']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths = [params.out_gro, params.out_cpt], raise_error = True, out = output.finished)

rule equil_complex_04_npt:
    input:
        top = approach_path + "/{ligand_name}/input/complex/complex.top",
        ndx = approach_path + "/{ligand_name}/input/complex/index.ndx",
        finished = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/03_npt/03_npt.finished",
        mdp = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/04_npt/04_npt.mdp",
    params:
        in_gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/03_npt/03_npt.gro",
        in_cpt = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/03_npt/03_npt.cpt",
        out_gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/04_npt/04_npt.gro",
        out_cpt = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/04_npt/04_npt.cpt",
        run_dir = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/04_npt",
    output:
        finished = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/04_npt/04_npt.finished",
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = params.in_gro,
            index = input.ndx,
            checkpoint = params.in_cpt,
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['complex']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths = [params.out_gro, params.out_cpt], raise_error = True, out = output.finished)

rule equil_complex_05_npt:
    input:
        top = approach_path + "/{ligand_name}/input/complex/complex.top",
        ndx = approach_path + "/{ligand_name}/input/complex/index.ndx",
        finished = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/04_npt/04_npt.finished",
        mdp = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/05_npt/05_npt.mdp",
    params:
        in_gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/04_npt/04_npt.gro",
        in_cpt = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/04_npt/04_npt.cpt",
        out_gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/05_npt/05_npt.gro",
        out_cpt = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/05_npt/05_npt.cpt",
        run_dir = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/05_npt",
    output:
        finished = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/05_npt/05_npt.finished",
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = params.in_gro,
            index = input.ndx,
            checkpoint = params.in_cpt,
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['complex']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths = [params.out_gro, params.out_cpt], raise_error = True, out = output.finished)

rule equil_complex_06_npt:
    input:
        top = approach_path + "/{ligand_name}/input/complex/complex.top",
        ndx = approach_path + "/{ligand_name}/input/complex/index.ndx",
        finished = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/05_npt/05_npt.finished",
        mdp = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/06_npt/06_npt.mdp",
    params:
        in_gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/05_npt/05_npt.gro",
        in_cpt = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/05_npt/05_npt.cpt",
        out_gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/06_npt/06_npt.gro",
        out_cpt = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/06_npt/06_npt.cpt",
        run_dir = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/06_npt"
    output:
        finished = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/06_npt/06_npt.finished",
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = params.in_gro,
            index = input.ndx,
            checkpoint = params.in_cpt,
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['complex']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths = [params.out_gro, params.out_cpt], raise_error = True, out = output.finished)

rule equil_complex_prod:
    input:
        top = approach_path + "/{ligand_name}/input/complex/complex.top",
        ndx = approach_path + "/{ligand_name}/input/complex/index.ndx",
        finished = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/06_npt/06_npt.finished",
        mdp = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.mdp",
    params:
        in_gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/06_npt/06_npt.gro",
        in_cpt = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/06_npt/06_npt.cpt",
        out_gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.gro",
        out_cpt = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.cpt",
        out_tpr = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.tpr",
        out_xtc = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.xtc",
        run_dir = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod",
    output:
        finished = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.finished",
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = params.in_gro,
            index = input.ndx,
            checkpoint = params.in_cpt,
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['complex']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths = [params.out_gro, params.out_cpt, params.out_tpr, params.out_xtc], raise_error = True, out = output.finished)