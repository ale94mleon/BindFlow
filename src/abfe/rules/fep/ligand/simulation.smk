from abfe.utils import tools

# Common to all the sub-workflows ligand/replica
approach_path = config["out_approach_path"]
threads = config['threads']
retries = config['retries']
load_dependencies = config['extra_directives']['dependencies']
mdrun_extra = config['extra_directives']['mdrun']

rule fep_ligand_00_min:
    input:
        top=approach_path + "/{ligand_name}/input/ligand/ligand.top",
        finished=approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/prod/prod.finished",
        mdp=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/00_min/00_min.mdp",
    params:
        in_gro=approach_path + "/{ligand_name}/{replica}/ligand/equil-mdsim/prod/prod.gro",
        run_dir=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/00_min",
    output:
        gro=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/00_min/00_min.gro"
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = params.in_gro,
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['ligand']
        )

rule fep_ligand_01_nvt:
    input:
        top=approach_path + "/{ligand_name}/input/ligand/ligand.top",
        gro=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/00_min/00_min.gro",
        mdp=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/01_nvt/01_nvt.mdp",
    params:
        out_gro=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/01_nvt/01_nvt.gro",
        out_cpt=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/01_nvt/01_nvt.cpt",
        run_dir=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/01_nvt",
    output:
        finished=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/01_nvt/01_nvt.finished",
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
        tools.paths_exist(paths = [params.out_gro, params.out_cpt], raise_error = True, out = output.finished)

rule fep_ligand_02_npt:
    input:
        top=approach_path + "/{ligand_name}/input/ligand/ligand.top",
        finished=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/01_nvt/01_nvt.finished",
        mdp=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/02_npt/02_npt.mdp",
    params:
        in_gro=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/01_nvt/01_nvt.gro",
        in_cpt=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/01_nvt/01_nvt.cpt",
        out_gro=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/02_npt/02_npt.gro",
        out_cpt=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/02_npt/02_npt.cpt",
        run_dir=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/02_npt",
    output:
        finished=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/02_npt/02_npt.finished",
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = params.in_gro,
            checkpoint = params.in_cpt,
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['ligand']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths = [params.out_gro, params.out_cpt], raise_error = True, out = output.finished)

rule fep_ligand_03_npt_norest:
    input:
        top=approach_path + "/{ligand_name}/input/ligand/ligand.top",
        finished=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/02_npt/02_npt.finished",
        mdp=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/03_npt_norest/03_npt_norest.mdp",
    params:
        in_gro=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/02_npt/02_npt.gro",
        in_cpt=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/02_npt/02_npt.cpt",
        out_gro=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/03_npt_norest/03_npt_norest.gro",
        out_cpt=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/03_npt_norest/03_npt_norest.cpt",
        run_dir=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/03_npt_norest",
    output:
        finished=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/03_npt_norest/03_npt_norest.finished",
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = params.in_gro,
            checkpoint = params.in_cpt,
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['ligand']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths = [params.out_gro, params.out_cpt], raise_error = True, out = output.finished)

rule fep_ligand_prod:
    input:
        top=approach_path + "/{ligand_name}/input/ligand/ligand.top",
        finished=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/03_npt_norest/03_npt_norest.finished",
        mdp=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/prod/prod.mdp",
    params:
        in_gro=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/03_npt_norest/03_npt_norest.gro",
        in_cpt=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/03_npt_norest/03_npt_norest.cpt",
        out_gro=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/prod/prod.gro",
        out_xvg=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/prod/prod.xvg",
        run_dir=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/prod",
    output:
        finished=approach_path + "/{ligand_name}/{replica}/ligand/fep/simulation/{state}/prod/prod.finished",
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = params.in_gro,
            checkpoint = params.in_cpt,
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['ligand']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths = [params.out_gro, params.out_xvg], raise_error = True, out = output.finished)