rule fep_complex_00_min:
    input:
        top=out_approach_path+"/{ligand_name}/{replica}/complex/fep/topology/complex_boresch.top",
        ndx=out_approach_path+"/{ligand_name}/{replica}/complex/fep/topology/index.ndx",
        gro=out_approach_path+"/{ligand_name}/{replica}/complex/equil-mdsim/boreschcalc/ClosestRestraintFrame.gro",
        mdp=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/00_min/00_min.mdp",
    params:
        run_dir=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/00_min/",
    output:
        gro=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/00_min/00_min.gro",
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

rule fep_complex_01_nvt:
    input:
        top=out_approach_path+"/{ligand_name}/{replica}/complex/fep/topology/complex_boresch.top",
        ndx=out_approach_path+"/{ligand_name}/{replica}/complex/fep/topology/index.ndx",
        gro=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/00_min/00_min.gro",
        mdp=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/01_nvt/01_nvt.mdp",
    params:
        out_gro=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/01_nvt/01_nvt.gro",
        out_cpt=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/01_nvt/01_nvt.cpt",
        run_dir=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/01_nvt",
    output:
        finished=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/01_nvt/01_nvt.finished",
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

rule fep_complex_02_npt:
    input:
        top=out_approach_path+"/{ligand_name}/{replica}/complex/fep/topology/complex_boresch.top",
        ndx=out_approach_path+"/{ligand_name}/{replica}/complex/fep/topology/index.ndx",
        finished=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/01_nvt/01_nvt.finished",
        mdp=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/02_npt/02_npt.mdp",
    params:
        in_gro=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/01_nvt/01_nvt.gro",
        in_cpt=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/01_nvt/01_nvt.cpt",
        out_gro=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/02_npt/02_npt.gro",
        out_cpt=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/02_npt/02_npt.cpt",
        run_dir=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/02_npt",
    output:
        finished=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/02_npt/02_npt.finished",
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

rule fep_complex_03_npt_norest:
    input:
        top=out_approach_path+"/{ligand_name}/{replica}/complex/fep/topology/complex_boresch.top",
        ndx=out_approach_path+"/{ligand_name}/{replica}/complex/fep/topology/index.ndx",
        finished=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/02_npt/02_npt.finished",
        mdp=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/03_npt_norest/03_npt_norest.mdp",
    params:
        in_gro=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/02_npt/02_npt.gro",
        in_cpt=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/02_npt/02_npt.cpt",
        out_gro=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/03_npt_norest/03_npt_norest.gro",
        out_cpt=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/03_npt_norest/03_npt_norest.cpt",
        run_dir=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/03_npt_norest",
    output:
        finished=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/03_npt_norest/03_npt_norest.finished",
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

rule fep_complex_prod:
    input:
        top=out_approach_path+"/{ligand_name}/{replica}/complex/fep/topology/complex_boresch.top",
        ndx=out_approach_path+"/{ligand_name}/{replica}/complex/fep/topology/index.ndx",
        finished=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/03_npt_norest/03_npt_norest.finished",
        mdp=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/prod/prod.mdp",
    params:
        in_gro=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/03_npt_norest/03_npt_norest.gro",
        in_cpt=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/03_npt_norest/03_npt_norest.cpt",
        out_gro=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/prod/prod.gro",
        out_xvg=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/prod/prod.xvg",
        run_dir=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/prod",
    output:
        finished=out_approach_path+"/{ligand_name}/{replica}/complex/fep/simulation/{state}/prod/prod.finished",
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
        tools.paths_exist(paths=[params.out_gro, params.out_xvg], raise_error=True, out=output.finished)