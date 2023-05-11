from abfe.utils.tools import gmx_runner

# Common to all the sub-workflows ligand/replica
run_path = config["run_path"]
input_path = config['input_data_path']
simulation_dir = run_path+"/ligand/fep/simulation"
threads = config['threads']
num_retries = config['num_retries']
load_dependencies = config['extra_directives']['dependencies']
mdrun_extra = config['extra_directives']['mdrun']

rule fep_run_ligand_emin:
    input:
        top=input_path+"/ligand/ligand.top",
        mdp=run_path+"/ligand/fep/simulation/{state}/emin/emin.mdp",
        gro=run_path+"/ligand/equil-mdsim/npt_equil2/npt_equil2.gro"
    params:
        run_dir=run_path+"/ligand/fep/simulation/{state}/emin",
    output:
        gro=run_path+"/ligand/fep/simulation/{state}/emin/emin.gro"
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
            **mdrun_extra
        )

rule fep_run_ligand_nvt_heat:
    input:
        top=input_path+"/ligand/ligand.top",
        mdp=run_path+"/ligand/fep/simulation/{state}/nvt/nvt.mdp",
        gro=run_path+"/ligand/fep/simulation/{state}/emin/emin.gro"
    params:
        run_dir=run_path+"/ligand/fep/simulation/{state}/nvt",
    output:
        gro=run_path+"/ligand/fep/simulation/{state}/nvt/nvt.gro",
        cpt=run_path+"/ligand/fep/simulation/{state}/nvt/nvt.cpt"
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
            **mdrun_extra
        )

rule fep_run_ligand_npt_eq1:
    input:
        top=input_path+"/ligand/ligand.top",
        mdp=run_path+"/ligand/fep/simulation/{state}/npt/npt.mdp",
        gro=run_path+"/ligand/fep/simulation/{state}/nvt/nvt.gro",
        cpt=run_path+"/ligand/fep/simulation/{state}/nvt/nvt.cpt"
    params:
        run_dir=run_path+"/ligand/fep/simulation/{state}/npt",
    output:
        gro=run_path+"/ligand/fep/simulation/{state}/npt/npt.gro",
        cpt=run_path+"/ligand/fep/simulation/{state}/npt/npt.cpt"
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
            **mdrun_extra
        )

rule fep_run_ligand_npt_eq2:
    input:
        top=input_path+"/ligand/ligand.top",
        mdp=run_path+"/ligand/fep/simulation/{state}/npt-norest/npt-norest.mdp",
        gro=run_path+"/ligand/fep/simulation/{state}/npt/npt.gro",
        cpt=run_path+"/ligand/fep/simulation/{state}/npt/npt.cpt"
    params:
        run_dir=run_path+"/ligand/fep/simulation/{state}/npt-norest",
    output:
        gro=run_path+"/ligand/fep/simulation/{state}/npt-norest/npt-norest.gro",
        cpt=run_path+"/ligand/fep/simulation/{state}/npt-norest/npt-norest.cpt"
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
            **mdrun_extra
        )

rule fep_run_ligand_prod:
    input:
        top=input_path+"/ligand/ligand.top",
        mdp=run_path+"/ligand/fep/simulation/{state}/prod/prod.mdp",
        gro=run_path+"/ligand/fep/simulation/{state}/npt-norest/npt-norest.gro",
        cpt=run_path+"/ligand/fep/simulation/{state}/npt-norest/npt-norest.cpt"
    params:
        run_dir=run_path+"/ligand/fep/simulation/{state}/prod",
    output:
        gro=run_path+"/ligand/fep/simulation/{state}/prod/prod.gro",
        xvg=run_path+"/ligand/fep/simulation/{state}/prod/prod.xvg"
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
            **mdrun_extra
        )
