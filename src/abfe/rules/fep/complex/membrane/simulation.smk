from abfe.utils.tools import gmx_runner


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
        mdp=run_path+"/complex/fep/simulation/{state}/00_min/00_min.mdp",
        gro=run_path+"/complex/equil-mdsim/boreschcalc/ClosestRestraintFrame.gro",
    params:
        run_dir=run_path+"/complex/fep/simulation/{state}/00_min/",
    output:
        gro=run_path+"/complex/fep/simulation/{state}/00_min/00_min.gro",
    threads: threads
    retries: retries
    run:
        gmx_runner(
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
        gro=run_path+"/complex/fep/simulation/{state}/00_min/00_min.gro"
    params:
        run_dir=run_path+"/complex/fep/simulation/{state}/01_nvt",
    output:
        gro=run_path+"/complex/fep/simulation/{state}/01_nvt/01_nvt.gro",
        cpt=run_path+"/complex/fep/simulation/{state}/01_nvt/01_nvt.cpt"
    threads: threads
    retries: retries
    run:
        gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = input.gro,
            index = input.ndx,
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['complex']
        )

rule fep_complex_02_npt:
    input:
        top=run_path+"/complex/fep/topology/complex_boresch.top",
        ndx=run_path+"/complex/fep/topology/index.ndx",
        mdp=run_path+"/complex/fep/simulation/{state}/02_npt/02_npt.mdp",
        gro=run_path+"/complex/fep/simulation/{state}/01_nvt/01_nvt.gro",
        cpt=run_path+"/complex/fep/simulation/{state}/01_nvt/01_nvt.cpt"
    params:
        run_dir=run_path+"/complex/fep/simulation/{state}/02_npt",
    output:
        gro=run_path+"/complex/fep/simulation/{state}/02_npt/02_npt.gro",
        cpt=run_path+"/complex/fep/simulation/{state}/02_npt/02_npt.cpt"
    threads: threads
    retries: retries
    run:
        gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = input.gro,
            index = input.ndx,
            checkpoint = input.cpt,
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['complex']
        )

rule fep_complex_03_npt_norest:
    input:
        top=run_path+"/complex/fep/topology/complex_boresch.top",
        ndx=run_path+"/complex/fep/topology/index.ndx",
        mdp=run_path+"/complex/fep/simulation/{state}/03_npt_norest/03_npt_norest.mdp",
        gro=run_path+"/complex/fep/simulation/{state}/02_npt/02_npt.gro",
        cpt=run_path+"/complex/fep/simulation/{state}/02_npt/02_npt.cpt"
    params:
        run_dir=run_path+"/complex/fep/simulation/{state}/03_npt_norest",
    output:
        gro=run_path+"/complex/fep/simulation/{state}/03_npt_norest/03_npt_norest.gro",
        cpt=run_path+"/complex/fep/simulation/{state}/03_npt_norest/03_npt_norest.cpt"
    threads: threads
    retries: retries
    run:
        gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = input.gro,
            index = input.ndx,
            checkpoint = input.cpt,
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['complex']
        )

rule fep_complex_prod:
    input:
        top=run_path+"/complex/fep/topology/complex_boresch.top",
        ndx=run_path+"/complex/fep/topology/index.ndx",
        mdp=run_path+"/complex/fep/simulation/{state}/prod/prod.mdp",
        gro=run_path+"/complex/fep/simulation/{state}/03_npt_norest/03_npt_norest.gro",
        cpt=run_path+"/complex/fep/simulation/{state}/03_npt_norest/03_npt_norest.cpt"
    params:
        run_dir=run_path+"/complex/fep/simulation/{state}/prod",
    output:
        gro=run_path+"/complex/fep/simulation/{state}/prod/prod.gro",
        xvg=run_path+"/complex/fep/simulation/{state}/prod/prod.xvg"
    threads: threads
    retries: retries
    run:
        gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = input.gro,
            index = input.ndx,
            checkpoint = input.cpt,
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['complex']
        )
