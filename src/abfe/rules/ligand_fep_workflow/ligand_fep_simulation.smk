from abfe.tools import gmx_runner

run_path = config["run_path"]
input_path = config['input_data_path']
simulation_dir = run_path+"/ligand/fep/simulation"

threads = config['threads']
num_retries = config['num_retries']
load_dependencies = config['job_extra_directives']
# TODO when is prepeared the configuration should be the dict with the key:values for mdrun or some dictionary 
# with the mdrun kwargs or an empty dict
mdrun_extra = config['mdrun_extra_directives']

rule fep_run_ligand_emin:
    input:
        top=input_path+"/ligand/ligand.top",
        mdp=simulation_dir+"/{state}/emin/emin.mdp",
        gro=run_path+"/ligand/equil-mdsim/npt_equil2/npt_equil2.gro"
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

rule fep_run_ligand_nvt_heat:
    input:
        top=input_path+"/ligand/ligand.top",
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
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            cpi = True,
            **mdrun_extra
        )

rule fep_run_ligand_npt_eq1:
    input:
        top=input_path+"/ligand/ligand.top",
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

rule fep_run_ligand_npt_eq2:
    input:
        top=input_path+"/ligand/ligand.top",
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
    shell:
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

rule fep_run_ligand_prod:
    input:
        top=input_path+"/ligand/ligand.top",
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
