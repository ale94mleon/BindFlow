from abfe.utils.tools import gmx_runner


# Common to all the sub-workflows ligand/replica
run_path = config["run_path"]
threads = config['threads']
num_retries = config['num_retries']
load_dependencies = config['extra_directives']['dependencies']
mdrun_extra = config['extra_directives']['mdrun']


rule fep_complex_emin:
    input:
        top=run_path+"/complex/fep/fep-topology/complex_boresch.top",
        mdp=run_path+"/complex/fep/simulation/{state}/emin/emin.mdp",
        gro=run_path+"/complex/equil-mdsim/boreschcalc/ClosestRestraintFrame.gro",
        
    params:
        run_dir=run_path+"/complex/fep/simulation/{state}/emin/",
    output:
        gro=run_path+"/complex/fep/simulation/{state}/emin/emin.gro",
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
            **mdrun_extra['complex']
        )

rule fep_complex_nvt:
    input:
        top=run_path+"/complex/fep/fep-topology/complex_boresch.top",
        mdp=run_path+"/complex/fep/simulation/{state}/nvt/nvt.mdp",
        gro=run_path+"/complex/fep/simulation/{state}/emin/emin.gro"
    params:
        run_dir=run_path+"/complex/fep/simulation/{state}/nvt",
    output:
        gro=run_path+"/complex/fep/simulation/{state}/nvt/nvt.gro",
        cpt=run_path+"/complex/fep/simulation/{state}/nvt/nvt.cpt"
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
            **mdrun_extra['complex']
        )

rule fep_complex_npt:
    input:
        top=run_path+"/complex/fep/fep-topology/complex_boresch.top",
        mdp=run_path+"/complex/fep/simulation/{state}/npt/npt.mdp",
        gro=run_path+"/complex/fep/simulation/{state}/nvt/nvt.gro",
        cpt=run_path+"/complex/fep/simulation/{state}/nvt/nvt.cpt"
    params:
        run_dir=run_path+"/complex/fep/simulation/{state}/npt",
    output:
        gro=run_path+"/complex/fep/simulation/{state}/npt/npt.gro",
        cpt=run_path+"/complex/fep/simulation/{state}/npt/npt.cpt"
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
            **mdrun_extra['complex']
        )

rule fep_complex_npt_norest:
    input:
        top=run_path+"/complex/fep/fep-topology/complex_boresch.top",
        mdp=run_path+"/complex/fep/simulation/{state}/npt_norest/npt_norest.mdp",
        gro=run_path+"/complex/fep/simulation/{state}/npt/npt.gro",
        cpt=run_path+"/complex/fep/simulation/{state}/npt/npt.cpt"
    params:
        run_dir=run_path+"/complex/fep/simulation/{state}/npt-norest",
    output:
        gro=run_path+"/complex/fep/simulation/{state}/npt_norest/npt_norest.gro",
        cpt=run_path+"/complex/fep/simulation/{state}/npt_norest/npt_norest.cpt"
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
            **mdrun_extra['complex']
        )

rule fep_complex_prod:
    input:
        top=run_path+"/complex/fep/fep-topology/complex_boresch.top",
        mdp=run_path+"/complex/fep/simulation/{state}/prod/prod.mdp",
        gro=run_path+"/complex/fep/simulation/{state}/npt_norest/npt_norest.gro",
        cpt=run_path+"/complex/fep/simulation/{state}/npt_norest/npt_norest.cpt"
    params:
        run_dir=run_path+"/complex/fep/simulation/{state}/prod",
    output:
        gro=run_path+"/complex/fep/simulation/{state}/prod/prod.gro",
        xvg=run_path+"/complex/fep/simulation/{state}/prod/prod.xvg"
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
            **mdrun_extra['complex']
        )
