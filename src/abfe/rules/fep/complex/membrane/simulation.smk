from abfe.utils import tools


# Common to all the sub-workflows ligand/replica
run_path = config["run_path"]
threads = config['threads']
retries = config['retries']
load_dependencies = config['extra_directives']['dependencies']
mdrun_extra = config['extra_directives']['mdrun']


rule fep_complex_00_min:
    input:
        top=run_path+"/complex/fep/topology/complex_boresch.top",
        ndx=run_path+"/complex/fep/topology/index.ndx",
        gro=run_path+"/complex/equil-mdsim/boreschcalc/ClosestRestraintFrame.gro",
        mdp=run_path+"/complex/fep/simulation/{state}/00_min/00_min.mdp", 
    params:
        run_dir=run_path+"/complex/fep/simulation/{state}/00_min/",
    output:
        gro=run_path+"/complex/fep/simulation/{state}/00_min/00_min.gro",
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

rule fep_complex_01_nvt:
    input:
        top=run_path+"/complex/fep/topology/complex_boresch.top",
        ndx=run_path+"/complex/fep/topology/index.ndx",
        mdp=run_path+"/complex/fep/simulation/{state}/01_nvt/01_nvt.mdp",
        gro=run_path+"/complex/fep/simulation/{state}/00_min/00_min.gro",
    params:
        out_gro=run_path+"/complex/fep/simulation/{state}/01_nvt/01_nvt.gro",
        out_cpt=run_path+"/complex/fep/simulation/{state}/01_nvt/01_nvt.cpt",
        run_dir=run_path+"/complex/fep/simulation/{state}/01_nvt",
    output:
        finished=run_path+"/complex/fep/simulation/{state}/01_nvt/01_nvt.finished",
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

rule fep_complex_02_npt:
    input:
        top=run_path+"/complex/fep/topology/complex_boresch.top",
        ndx=run_path+"/complex/fep/topology/index.ndx",
        finished=run_path+"/complex/fep/simulation/{state}/01_nvt/01_nvt.finished",
        mdp=run_path+"/complex/fep/simulation/{state}/02_npt/02_npt.mdp",
    params:
        in_gro=run_path+"/complex/fep/simulation/{state}/01_nvt/01_nvt.gro",
        in_cpt=run_path+"/complex/fep/simulation/{state}/01_nvt/01_nvt.cpt",
        out_gro=run_path+"/complex/fep/simulation/{state}/02_npt/02_npt.gro",
        out_cpt=run_path+"/complex/fep/simulation/{state}/02_npt/02_npt.cpt",
        run_dir=run_path+"/complex/fep/simulation/{state}/02_npt",
    output:
        finished=run_path+"/complex/fep/simulation/{state}/02_npt/02_npt.finished",
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

rule fep_complex_03_npt_norest:
    input:
        top=run_path+"/complex/fep/topology/complex_boresch.top",
        ndx=run_path+"/complex/fep/topology/index.ndx",
        finished=run_path+"/complex/fep/simulation/{state}/02_npt/02_npt.finished",
        mdp=run_path+"/complex/fep/simulation/{state}/03_npt_norest/03_npt_norest.mdp",
    params:
        in_gro=run_path+"/complex/fep/simulation/{state}/02_npt/02_npt.gro",
        in_cpt=run_path+"/complex/fep/simulation/{state}/02_npt/02_npt.cpt",
        out_gro=run_path+"/complex/fep/simulation/{state}/03_npt_norest/03_npt_norest.gro",
        out_cpt=run_path+"/complex/fep/simulation/{state}/03_npt_norest/03_npt_norest.cpt",
        run_dir=run_path+"/complex/fep/simulation/{state}/03_npt_norest",
    output:
        finished=run_path+"/complex/fep/simulation/{state}/03_npt_norest/03_npt_norest.finished",
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

rule fep_complex_prod:
    input:
        top=run_path+"/complex/fep/topology/complex_boresch.top",
        ndx=run_path+"/complex/fep/topology/index.ndx",
        finished=run_path+"/complex/fep/simulation/{state}/03_npt_norest/03_npt_norest.finished",
        mdp=run_path+"/complex/fep/simulation/{state}/prod/prod.mdp",  
    params:
        in_gro=run_path+"/complex/fep/simulation/{state}/03_npt_norest/03_npt_norest.gro",
        in_cpt=run_path+"/complex/fep/simulation/{state}/03_npt_norest/03_npt_norest.cpt",
        out_gro=run_path+"/complex/fep/simulation/{state}/prod/prod.gro",
        out_xvg=run_path+"/complex/fep/simulation/{state}/prod/prod.xvg",
        run_dir=run_path+"/complex/fep/simulation/{state}/prod",
    output:
        finished=run_path+"/complex/fep/simulation/{state}/prod/prod.finished",
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
        tools.paths_exist(paths = [params.out_gro, params.out_xvg], raise_error = True, out = output.finished)