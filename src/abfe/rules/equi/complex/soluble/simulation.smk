from abfe.utils import tools
from abfe.preparation import boresch
from abfe.mdp import mdp

# Common to all the sub-workflows ligand/replica
out_approach_path = config["out_approach_path"]
threads = config['threads']
retries = config['retries']
load_dependencies = config['extra_directives']['dependencies']
mdrun_extra = config['extra_directives']['mdrun']
rule equil_complex_00_min:
    input:
        top = approach_path + "/{ligand_name}/input/complex/complex.top",
        gro = approach_path + "/{ligand_name}/input/complex/complex.gro",
        mdp = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/00_min/00_min.mdp"
    params:
        run_dir = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/00_min"
    output:
        gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/00_min/00_min.gro"
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

rule equil_complex_01_nvt:
    input:
        top = approach_path + "/{ligand_name}/input/complex/complex.top",
        gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/00_min/00_min.gro",
        mdp = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt/01_nvt.mdp",
    params:
        out_gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt/01_nvt.gro",
        out_cpt = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt/01_nvt.cpt",
        run_dir = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt",
    output:
        finished = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt/01_nvt.finished",
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
        tools.paths_exist(paths = [params.out_gro, params.out_cpt], raise_error = True, out = output.finished)

rule equil_complex_02_nvt:
    input:
        top = approach_path + "/{ligand_name}/input/complex/complex.top",
        finished = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt/01_nvt.finished",
        mdp = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/02_nvt/02_nvt.mdp",
    params:
        in_gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt/01_nvt.gro",
        in_cpt = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt/01_nvt.cpt",
        out_gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/02_nvt/02_nvt.gro",
        out_cpt = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/02_nvt/02_nvt.cpt",
        run_dir = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/02_nvt",
    output:
        finished = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/02_nvt/02_nvt.finished",
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
            **mdrun_extra['complex']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths = [params.out_gro, params.out_cpt], raise_error = True, out = output.finished)

rule equil_complex_03_npt:
    input:
        top = approach_path + "/{ligand_name}/input/complex/complex.top",
        finished = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/02_nvt/02_nvt.finished",
        mdp = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/03_npt/03_npt.mdp",
    params:
        in_gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/02_nvt/02_nvt.gro",
        in_cpt = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/02_nvt/02_nvt.cpt",
        out_gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/03_npt/03_npt.gro",
        out_cpt = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/03_npt/03_npt.cpt",
        run_dir = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/03_npt",
    output:
        finished=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/03_npt/03_npt.finished",
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
            **mdrun_extra['complex']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths = [params.out_gro, params.out_cpt], raise_error = True, out = output.finished)

rule equil_complex_04_npt:
    input:
        top = approach_path + "/{ligand_name}/input/complex/complex.top",
        finished = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/03_npt/03_npt.finished",
        mdp = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/04_npt/04_npt.mdp",
    params:
        in_gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/03_npt/03_npt.gro",
        in_cpt = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/03_npt/03_npt.cpt",
        out_gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/04_npt/04_npt.gro",
        out_cpt = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/04_npt/04_npt.cpt",
        run_dir = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/04_npt",
    output:
        finished=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/04_npt/04_npt.finished",
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
            **mdrun_extra['complex']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths = [params.out_gro, params.out_cpt], raise_error = True, out = output.finished)

rule equil_complex_prod:
    input:
        top = approach_path + "/{ligand_name}/input/complex/complex.top",
        finished = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/04_npt/04_npt.finished",
        mdp = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.mdp",
    params:
        in_gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/04_npt/04_npt.gro",
        in_cpt = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/04_npt/04_npt.cpt",
        out_gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.gro",
        out_cpt = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.cpt",
        out_tpr = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.tpr",
        out_xtc = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.xtc",
        run_dir = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod",
    output:
        finished = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.finished",
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
            **mdrun_extra['complex']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths = [params.out_gro, params.out_cpt, params.out_tpr, params.out_xtc], raise_error = True, out = output.finished)

rule equil_complex_get_boresch_restraints:
    input:
        finished = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.finished",
        mdp = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.mdp",
    params:
        in_tpr = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.tpr",
        in_xtc = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.xtc",
        run_dir = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/boreschcalc",
    output:
        gro = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/boreschcalc/ClosestRestraintFrame.gro",
        top = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/boreschcalc/BoreschRestraint.top",
        boresch_dG_off = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/boreschcalc/dG_off.dat"
    run:
        # Fix trajectory.
        tools.makedirs(params.run_dir)
        tools.run(f"export GMX_MAXBACKUP=-1; echo 'System' | gmx trjconv -s {params.in_tpr} -f {params.in_xtc} -o {params.run_dir}/whole.xtc -pbc whole")
        tools.run(f"export GMX_MAXBACKUP=-1; echo 'System' | gmx trjconv -s {params.in_tpr} -f {params.run_dir}/whole.xtc -o {params.run_dir}/nojump.xtc -pbc nojump")
        tools.run(f"export GMX_MAXBACKUP=-1; echo 'Protein System' | gmx trjconv -s {params.in_tpr} -f {params.run_dir}/nojump.xtc -o {params.run_dir}/prod_center.xtc -pbc mol -center -ur compact")
        # Clean
        tools.run(f"rm {params.run_dir}/whole.xtc {params.run_dir}/nojump.xtc")

        # Getting Borech restraints
        mdp_params = mdp.MDP().from_file(input.mdp).parameters
        if 'ref-t' in mdp_params:
            temperature = float(mdp_params['ref-t'].split()[0])
        elif 'ref_t' in mdp_params:
            temperature = float(mdp_params['ref_t'].split()[0])
        boresch.gen_restraint(
            topology = params.in_tpr,
            trajectory = f"{params.run_dir}/prod_center.xtc",
            outpath = params.run_dir,
            temperature = temperature
        )
        # Clean
        tools.run(f"rm {params.run_dir}/prod_center.xtc")