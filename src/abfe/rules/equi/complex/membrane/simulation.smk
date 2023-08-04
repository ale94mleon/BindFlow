from abfe.utils import tools
from abfe.preparation import boresch
from abfe.mdp import mdp

# Common to all the sub-workflows ligand/replica
approach_path = config["out_approach_path"]
threads = config['threads']
retries = config['retries']
load_dependencies = config['extra_directives']['dependencies']
mdrun_extra = config['extra_directives']['mdrun']

rule equil_complex_00_min:
    input:
        top=approach_path + "/{ligand_name}/input/complex/complex.top",
        ndx=approach_path + "/{ligand_name}/input/complex/index.ndx",
        gro=approach_path + "/{ligand_name}/input/complex/complex.gro",
        mdp=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/00_min/00_min.mdp"
    params:
        run_dir=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/00_min"
    output:
        gro=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/00_min/00_min.gro"
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

rule equil_complex_01_nvt:
    input:
        top=approach_path + "/{ligand_name}/input/complex/complex.top",
        ndx=approach_path + "/{ligand_name}/input/complex/index.ndx",
        gro=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/00_min/00_min.gro",
        mdp=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt/01_nvt.mdp",
    params:
        output={
            'gro':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt/01_nvt.gro",
            'cpt':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt/01_nvt.cpt"
        },
        run_dir=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt",
    output:
        finished=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt/01_nvt.finished",
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
        tools.paths_exist(paths = params.output.values(), raise_error = True, out = output.finished)

rule equil_complex_02_npt:
    input:
        top=approach_path + "/{ligand_name}/input/complex/complex.top",
        ndx=approach_path + "/{ligand_name}/input/complex/index.ndx",
        finished=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt/01_nvt.finished",
        mdp=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/02_npt/02_npt.mdp",
    params:
        input={
            'gro':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt/01_nvt.gro",
            'cpt':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/01_nvt/01_nvt.cpt"
        },
        output={
            'gro':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/02_npt/02_npt.gro",
            'cpt':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/02_npt/02_npt.cpt"
        },
        run_dir=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/02_npt",
    output:
        finished=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/02_npt/02_npt.finished",
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = params.input['gro'],
            index = input.ndx,
            checkpoint = params.input['cpt'],
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['complex']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths = params.output.values(), raise_error = True, out = output.finished)

rule equil_complex_03_npt:
    input:
        top=approach_path + "/{ligand_name}/input/complex/complex.top",
        ndx=approach_path + "/{ligand_name}/input/complex/index.ndx",
        finished=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/02_npt/02_npt.finished",
        mdp=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/03_npt/03_npt.mdp",
    params:
        input={
            'gro':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/02_npt/02_npt.gro",
            'cpt':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/02_npt/02_npt.cpt"
        },
        output={
            'gro':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/03_npt/03_npt.gro",
            'cpt':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/03_npt/03_npt.cpt"
        },
        run_dir=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/03_npt",
    output:
        finished=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/03_npt/03_npt.finished",
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = params.input['gro'],
            index = input.ndx,
            checkpoint = params.input['cpt'],
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['complex']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths = params.output.values(), raise_error = True, out = output.finished)

rule equil_complex_04_npt:
    input:
        top=approach_path + "/{ligand_name}/input/complex/complex.top",
        ndx=approach_path + "/{ligand_name}/input/complex/index.ndx",
        finished=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/03_npt/03_npt.finished",
        mdp=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/04_npt/04_npt.mdp",
    params:
        input={
            'gro':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/03_npt/03_npt.gro",
            'cpt':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/03_npt/03_npt.cpt"
        },
        output={
            'gro':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/04_npt/04_npt.gro",
            'cpt':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/04_npt/04_npt.cpt"
        },
        run_dir=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/04_npt"
    output:
        finished=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/04_npt/04_npt.finished",
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = params.input['gro'],
            index = input.ndx,
            checkpoint = params.input['cpt'],
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['complex']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths = params.output.values(), raise_error = True, out = output.finished)

rule equil_complex_05_npt:
    input:
        top=approach_path + "/{ligand_name}/input/complex/complex.top",
        ndx=approach_path + "/{ligand_name}/input/complex/index.ndx",
        finished=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/04_npt/04_npt.finished",
        mdp=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/05_npt/05_npt.mdp",
    params:
        input={
            'gro':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/04_npt/04_npt.gro",
            'cpt':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/04_npt/04_npt.cpt"
        },
        output={
            'gro':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/05_npt/05_npt.gro",
            'cpt':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/05_npt/05_npt.cpt"
        },
        run_dir=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/05_npt"
    output:
        finished=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/05_npt/05_npt.finished",
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = params.input['gro'],
            index = input.ndx,
            checkpoint = params.input['cpt'],
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['complex']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths = params.output.values(), raise_error = True, out = output.finished)

rule equil_complex_06_npt:
    input:
        top=approach_path + "/{ligand_name}/input/complex/complex.top",
        ndx=approach_path + "/{ligand_name}/input/complex/index.ndx",
        finished=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/05_npt/05_npt.finished",
        mdp=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/06_npt/06_npt.mdp",
    params:
        input={
            'gro':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/05_npt/05_npt.gro",
            'cpt':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/05_npt/05_npt.cpt"
        },
        output={
            'gro':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/06_npt/06_npt.gro",
            'cpt':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/06_npt/06_npt.cpt"
        },
        run_dir=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/06_npt"
    output:
        finished=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/06_npt/06_npt.finished",
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = params.input['gro'],
            index = input.ndx,
            checkpoint = params.input['cpt'],
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['complex']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths = params.output.values(), raise_error = True, out = output.finished)

rule equil_complex_prod:
    input:
        top=approach_path + "/{ligand_name}/input/complex/complex.top",
        ndx=approach_path + "/{ligand_name}/input/complex/index.ndx",
        finished=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/06_npt/06_npt.finished",
        mdp=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.mdp",
    params:
        input={
            'gro':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/06_npt/06_npt.gro",
            'cpt':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/06_npt/06_npt.cpt"
        },
        output={
            'gro':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.gro",
            'cpt':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.cpt",
            'tpr':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.tpr",
            'xtc':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.xtc",
        },
        run_dir=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod",
    output:
        finished=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.finished",
    threads: threads
    retries: retries
    run:
        tools.gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = params.input['gro'],
            index = input.ndx,
            checkpoint = params.input['cpt'],
            nthreads = threads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            **mdrun_extra['complex']
        )
        # Allow proper GROMACS continuation
        tools.paths_exist(paths = params.output.values(), raise_error = True, out = output.finished)

rule equil_run_complex_trjconv:
    input:
        finished=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.finished",
    params:
        input={
            'tpr':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.tpr",
            'xtc':approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.xtc"
        },
        run_dir=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/boreschcalc/"
    output:
        xtc=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/boreschcalc/prod_center.xtc"
    run:
        tools.makedirs(params.run_dir)
        tools.run(f"echo 'System' | gmx trjconv -s {params.input['tpr']} -f {params.input['xtc']} -o {params.run_dir}/whole.xtc -pbc whole")
        tools.run(f"echo 'System' | gmx trjconv -s {params.input['tpr']} -f {params.run_dir}/whole.xtc -o {params.run_dir}/nojump.xtc -pbc nojump")
        tools.run(f"echo 'Protein System' | gmx trjconv -s {params.input['tpr']} -f {params.run_dir}/nojump.xtc -o {output.xtc} -pbc mol -center -ur compact")
        tools.run(f"rm {params.run_dir}/whole.xtc {params.run_dir}/nojump.xtc")

rule equil_run_complex_get_boresch_restraints:
    input:
        finished=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.finished",
        xtc=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/boreschcalc/prod_center.xtc",
        mdp=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.mdp",
    params:
        tpr=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.tpr",
        run_dir=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/boreschcalc/",
    output:
        gro=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/boreschcalc/ClosestRestraintFrame.gro",
        top=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/boreschcalc/BoreschRestraint.top",
        boresch_dG_off=approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/boreschcalc/dG_off.dat"
    run:
        mdp_params = mdp.MDP().from_file(input.mdp).parameters
        if 'ref-t' in mdp_params:
            temperature = float(mdp_params['ref-t'].split()[0])
        elif 'ref_t' in mdp_params:
            temperature = float(mdp_params['ref_t'].split()[0])
        boresch.gen_restraint(
            topology = params.tpr,
            trajectory = input.xtc,
            outpath = params.run_dir,
            temperature = temperature
        )
