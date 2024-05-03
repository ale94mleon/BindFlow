from bindflow.utils import tools


# Common to all the sub-workflows ligand/replica
approach_path = config["out_approach_path"]
threads = config['threads']
retries = config['retries']
load_dependencies = config['extra_directives']['dependencies']
mdrun_extra = config['extra_directives']['mdrun']


rule mmxbsa_sample_prod:
    input:
        top = approach_path + "/{ligand_name}/input/complex/complex.top",
        mdp = approach_path + "/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/prod.mdp",
        gro = approach_path + "/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/init.gro",
    params:
        out_gro = approach_path + "/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/prod.gro",
        out_cpt = approach_path + "/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/prod.cpt",
        out_tpr = approach_path + "/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/prod.tpr",
        out_xtc = approach_path + "/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/prod.xtc",
        run_dir = approach_path + "/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/",
    output:
        finished =  approach_path + "/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/prod.finished",
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
            **mdrun_extra['complex']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths = [params.out_gro, params.out_cpt, params.out_tpr, params.out_xtc], raise_error = True, out = output.finished)
