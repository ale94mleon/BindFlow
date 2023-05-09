from abfe.utils.tools import gmx_runner

run_path = config["run_path"]
simulation_dir = run_path+"/complex/fep/simulation"

threads = config['threads']
num_retries = config['num_retries']
load_dependencies = config['extra_directives']['dependencies']
mdrun_extra = config['extra_directives']['mdrun']


rule fep_run_complex_emin:
    input:
        top=run_path+"/complex/fep/fep-topology/complex_boresch.top",
        mdp=simulation_dir+"/{state}/emin/emin.mdp",
        gro=run_path+"/complex/equil-mdsim/boreschcalc/ClosestRestraintFrame.gro",
        
    params:
        run_dir=simulation_dir+"/{state}/emin",
    output:
        gro=simulation_dir+"/{state}/emin/emin.gro"
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
            cpi = True,
            **mdrun_extra
        )

rule fep_run_complex_nvt_heat:
    input:
        top=run_path+"/complex/fep/fep-topology/complex_boresch.top",
        mdp=simulation_dir+"/{state}/nvt/nvt.mdp",
        gro=simulation_dir+"/{state}/emin/emin.gro"
    params:
        run_dir=simulation_dir+"/{state}/nvt",
    output:
        gro=simulation_dir+"/{state}/nvt/nvt.gro",
        cpt=simulation_dir+"/{state}/nvt/nvt.cpt"
    threads: threads
    retries: num_retries
    run:
        gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = input.gro,
            nthreads = threadss,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            cpi = True,
            **mdrun_extra
        )

rule fep_run_complex_npt_eq1:
    input:
        top=run_path+"/complex/fep/fep-topology/complex_boresch.top",
        mdp=simulation_dir+"/{state}/npt/npt.mdp",
        gro=simulation_dir+"/{state}/nvt/nvt.gro",
        cpt=simulation_dir+"/{state}/nvt/nvt.cpt"
    params:
        run_dir=simulation_dir+"/{state}/npt",
    output:
        gro=simulation_dir+"/{state}/npt/npt.gro",
        cpt=simulation_dir+"/{state}/npt/npt.cpt"
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
            cpi = True,
            **mdrun_extra
        )

rule fep_run_complex_npt_eq2:
    input:
        top=run_path+"/complex/fep/fep-topology/complex_boresch.top",
        mdp=simulation_dir+"/{state}/npt-norest/npt-norest.mdp",
        gro=simulation_dir+"/{state}/npt/npt.gro",
        cpt=simulation_dir+"/{state}/npt/npt.cpt"
    params:
        run_dir=simulation_dir+"/{state}/npt-norest",
    output:
        gro=simulation_dir+"/{state}/npt-norest/npt-norest.gro",
        cpt=simulation_dir+"/{state}/npt-norest/npt-norest.cpt"
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
            cpi = True,
            **mdrun_extra
        )

rule fep_run_complex_prod:
    input:
        top=run_path+"/complex/fep/fep-topology/complex_boresch.top",
        mdp=simulation_dir+"/{state}/prod/prod.mdp",
        gro=simulation_dir+"/{state}/npt-norest/npt-norest.gro",
        cpt=simulation_dir+"/{state}/npt-norest/npt-norest.cpt"
    params:
        run_dir=simulation_dir+"/{state}/prod",
    output:
        gro=simulation_dir+"/{state}/prod/prod.gro",
        xvg=simulation_dir+"/{state}/prod/prod.xvg"
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
            cpi = True,
            **mdrun_extra
        )
