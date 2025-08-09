rule equil_complex_00_min:
    input:
        top=out_approach_path+"/{ligand_name}/input/complex/complex.top",
        ndx=out_approach_path+"/{ligand_name}/input/complex/index.ndx",
        gro=out_approach_path+"/{ligand_name}/input/complex/complex.gro",
        mdp=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/00_min/00_min.mdp"
    params:
        run_dir=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/00_min"
    output:
        gro=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/00_min/00_min.gro"
    threads: threads
    retries: retries
    run:
        update_mdrun_extra = mdrun_extra["complex"].copy()
        for invalid_flag in ["update", "bonded"]:
            if invalid_flag in update_mdrun_extra:
                del update_mdrun_extra[invalid_flag]
        tools.gmx_runner(
            mdp=input.mdp,
            topology=input.top,
            structure=input.gro,
            index=input.ndx,
            nthreads=threads,
            load_dependencies=load_dependencies,
            run_dir=params.run_dir,
            **update_mdrun_extra
        )

rule equil_complex_01_nvt:
    input:
        top=out_approach_path+"/{ligand_name}/input/complex/complex.top",
        ndx=out_approach_path+"/{ligand_name}/input/complex/index.ndx",
        gro=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/00_min/00_min.gro",
        mdp=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt/01_nvt.mdp",
    params:
        out_gro=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt/01_nvt.gro",
        out_cpt=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt/01_nvt.cpt",
        run_dir=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt",
    output:
        finished=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt/01_nvt.finished",
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp=input.mdp,
            topology=input.top,
            structure=input.gro,
            index=input.ndx,
            nthreads=threads,
            load_dependencies=load_dependencies,
            run_dir=params.run_dir,
            **mdrun_extra['complex']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths=[params.out_gro, params.out_cpt], raise_error=True, out=output.finished)

rule equil_complex_02_nvt:
    input:
        top=out_approach_path+"/{ligand_name}/input/complex/complex.top",
        ndx=out_approach_path+"/{ligand_name}/input/complex/index.ndx",
        finished=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt/01_nvt.finished",
        mdp=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/02_nvt/02_nvt.mdp",
    params:
        in_gro=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt/01_nvt.gro",
        in_cpt=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt/01_nvt.cpt",
        out_gro=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/02_nvt/02_nvt.gro",
        out_cpt=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/02_nvt/02_nvt.cpt",
        run_dir=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/02_nvt",
    output:
        finished=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/02_nvt/02_nvt.finished",
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp=input.mdp,
            topology=input.top,
            structure=params.in_gro,
            index=input.ndx,
            checkpoint=params.in_cpt,
            nthreads=threads,
            load_dependencies=load_dependencies,
            run_dir=params.run_dir,
            **mdrun_extra['complex']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths=[params.out_gro, params.out_cpt], raise_error=True, out=output.finished)

rule equil_complex_03_npt:
    input:
        top=out_approach_path+"/{ligand_name}/input/complex/complex.top",
        ndx=out_approach_path+"/{ligand_name}/input/complex/index.ndx",
        finished=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/02_nvt/02_nvt.finished",
        mdp=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/03_npt/03_npt.mdp",
    params:
        in_gro=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/02_nvt/02_nvt.gro",
        in_cpt=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/02_nvt/02_nvt.cpt",
        out_gro=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/03_npt/03_npt.gro",
        out_cpt=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/03_npt/03_npt.cpt",
        run_dir=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/03_npt",
    output:
        finished=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/03_npt/03_npt.finished",
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp=input.mdp,
            topology=input.top,
            structure=params.in_gro,
            index=input.ndx,
            checkpoint=params.in_cpt,
            nthreads=threads,
            load_dependencies=load_dependencies,
            run_dir=params.run_dir,
            **mdrun_extra['complex']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths=[params.out_gro, params.out_cpt], raise_error=True, out=output.finished)

rule equil_complex_04_npt:
    input:
        top=out_approach_path+"/{ligand_name}/input/complex/complex.top",
        ndx=out_approach_path+"/{ligand_name}/input/complex/index.ndx",
        finished=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/03_npt/03_npt.finished",
        mdp=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/04_npt/04_npt.mdp",
    params:
        in_gro=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/03_npt/03_npt.gro",
        in_cpt=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/03_npt/03_npt.cpt",
        out_gro=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/04_npt/04_npt.gro",
        out_cpt=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/04_npt/04_npt.cpt",
        run_dir=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/04_npt",
    output:
        finished=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/04_npt/04_npt.finished",
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp=input.mdp,
            topology=input.top,
            structure=params.in_gro,
            index=input.ndx,
            checkpoint=params.in_cpt,
            nthreads=threads,
            load_dependencies=load_dependencies,
            run_dir=params.run_dir,
            **mdrun_extra['complex']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths=[params.out_gro, params.out_cpt], raise_error=True, out=output.finished)

rule equil_complex_05_npt:
    input:
        top=out_approach_path+"/{ligand_name}/input/complex/complex.top",
        ndx=out_approach_path+"/{ligand_name}/input/complex/index.ndx",
        finished=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/04_npt/04_npt.finished",
        mdp=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/05_npt/05_npt.mdp",
    params:
        in_gro=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/04_npt/04_npt.gro",
        in_cpt=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/04_npt/04_npt.cpt",
        out_gro=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/05_npt/05_npt.gro",
        out_cpt=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/05_npt/05_npt.cpt",
        run_dir=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/05_npt",
    output:
        finished=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/05_npt/05_npt.finished",
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp=input.mdp,
            topology=input.top,
            structure=params.in_gro,
            index=input.ndx,
            checkpoint=params.in_cpt,
            nthreads=threads,
            load_dependencies=load_dependencies,
            run_dir=params.run_dir,
            **mdrun_extra['complex']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths=[params.out_gro, params.out_cpt], raise_error=True, out=output.finished)

rule equil_complex_06_npt:
    input:
        top=out_approach_path+"/{ligand_name}/input/complex/complex.top",
        ndx=out_approach_path+"/{ligand_name}/input/complex/index.ndx",
        finished=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/05_npt/05_npt.finished",
        mdp=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/06_npt/06_npt.mdp",
    params:
        in_gro=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/05_npt/05_npt.gro",
        in_cpt=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/05_npt/05_npt.cpt",
        out_gro=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/06_npt/06_npt.gro",
        out_cpt=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/06_npt/06_npt.cpt",
        run_dir=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/06_npt"
    output:
        finished=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/06_npt/06_npt.finished",
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp=input.mdp,
            topology=input.top,
            structure=params.in_gro,
            index=input.ndx,
            checkpoint=params.in_cpt,
            nthreads=threads,
            load_dependencies=load_dependencies,
            run_dir=params.run_dir,
            **mdrun_extra['complex']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths=[params.out_gro, params.out_cpt], raise_error=True, out=output.finished)

rule equil_complex_prod:
    input:
        top=out_approach_path+"/{ligand_name}/input/complex/complex.top",
        ndx=out_approach_path+"/{ligand_name}/input/complex/index.ndx",
        finished=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/06_npt/06_npt.finished",
        mdp=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.mdp",
    params:
        in_gro=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/06_npt/06_npt.gro",
        in_cpt=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/06_npt/06_npt.cpt",
        out_gro=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.gro",
        out_cpt=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.cpt",
        out_tpr=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.tpr",
        out_xtc=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.xtc",
        run_dir=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/prod",
    output:
        finished=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.finished",
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp=input.mdp,
            topology=input.top,
            structure=params.in_gro,
            index=input.ndx,
            checkpoint=params.in_cpt,
            nthreads=threads,
            load_dependencies=load_dependencies,
            run_dir=params.run_dir,
            **mdrun_extra['complex']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths = [params.out_gro, params.out_cpt, params.out_tpr, params.out_xtc], raise_error = True, out = output.finished)