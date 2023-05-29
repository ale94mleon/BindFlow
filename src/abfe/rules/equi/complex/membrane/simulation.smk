from abfe.utils.tools import gmx_runner
from abfe.scripts.preparation import boresch

# Common to all the sub-workflows ligand/replica
input_path = config['input_data_path']
run_path = config["run_path"]
threads = config['threads']
num_retries = config['num_retries']
load_dependencies = config['extra_directives']['dependencies']
mdrun_extra = config['extra_directives']['mdrun']

rule equil_complex_min0:
    input:
        top=input_path+"/complex/complex.top",
        ndx=input_path+"/complex/index.ndx",
        gro=input_path+"/complex/complex.gro",
        mdp=run_path+"/complex/equil-mdsim/00_min/00_min.mdp"
    params:
        run_dir=run_path+"/complex/equil-mdsim/00_min"
    output:
        gro=run_path+"/complex/equil-mdsim/00_min/00_min.gro"
    threads: threads
    retries: num_retries
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

rule equil_complex_nvt1:
    input:
        top=input_path+"/complex/complex.top",
        ndx=input_path+"/complex/index.ndx",
        gro=run_path+"/complex/equil-mdsim/00_min/00_min.gro",
        mdp=run_path+"/complex/equil-mdsim/01_nvt/01_nvt.mdp",
    params:
        run_dir=run_path+"/complex/equil-mdsim/01_nvt"
    output:
        gro=run_path+"/complex/equil-mdsim/01_nvt/01_nvt.gro",
        cpt=run_path+"/complex/equil-mdsim/01_nvt/01_nvt.cpt"
    threads: threads
    retries: num_retries
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

rule equil_complex_npt2:
    input:
        top=input_path+"/complex/complex.top",
        ndx=input_path+"/complex/index.ndx",
        gro=run_path+"/complex/equil-mdsim/01_nvt/01_nvt.gro",
        cpt=run_path+"/complex/equil-mdsim/01_nvt/01_nvt.cpt",
        mdp=run_path+"/complex/equil-mdsim/02_npt/02_npt.mdp",
    params:
        run_dir=run_path+"/complex/equil-mdsim/02_npt"
    output:
        gro=run_path+"/complex/equil-mdsim/02_npt/02_npt.gro",
        cpt=run_path+"/complex/equil-mdsim/02_npt/02_npt.cpt"
    threads: threads
    retries: num_retries
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

rule equil_complex_npt3:
    input:
        top=input_path+"/complex/complex.top",
        ndx=input_path+"/complex/index.ndx",
        gro=run_path+"/complex/equil-mdsim/02_npt/02_npt.gro",
        cpt=run_path+"/complex/equil-mdsim/02_npt/02_npt.cpt",
        mdp=run_path+"/complex/equil-mdsim/03_npt/03_npt.mdp",
    params:
        run_dir=run_path+"/complex/equil-mdsim/03_npt"
    output:
        gro=run_path+"/complex/equil-mdsim/03_npt/03_npt.gro",
        cpt=run_path+"/complex/equil-mdsim/03_npt/03_npt.cpt"
    threads: threads
    retries: num_retries
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

rule equil_complex_npt4:
    input:
        top=input_path+"/complex/complex.top",
        ndx=input_path+"/complex/index.ndx",
        gro=run_path+"/complex/equil-mdsim/03_npt/03_npt.gro",
        cpt=run_path+"/complex/equil-mdsim/03_npt/03_npt.cpt",
        mdp=run_path+"/complex/equil-mdsim/04_npt/04_npt.mdp",
    params:
        run_dir=run_path+"/complex/equil-mdsim/04_npt"
    output:
        gro=run_path+"/complex/equil-mdsim/04_npt/04_npt.gro",
        cpt=run_path+"/complex/equil-mdsim/04_npt/04_npt.cpt"
    threads: threads
    retries: num_retries
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

rule equil_complex_npt5:
    input:
        top=input_path+"/complex/complex.top",
        ndx=input_path+"/complex/index.ndx",
        gro=run_path+"/complex/equil-mdsim/04_npt/04_npt.gro",
        cpt=run_path+"/complex/equil-mdsim/04_npt/04_npt.cpt",
        mdp=run_path+"/complex/equil-mdsim/05_npt/05_npt.mdp",
    params:
        run_dir=run_path+"/complex/equil-mdsim/05_npt"
    output:
        gro=run_path+"/complex/equil-mdsim/05_npt/05_npt.gro",
        cpt=run_path+"/complex/equil-mdsim/05_npt/05_npt.cpt"
    threads: threads
    retries: num_retries
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
rule equil_complex_npt6:
    input:
        top=input_path+"/complex/complex.top",
        ndx=input_path+"/complex/index.ndx",
        gro=run_path+"/complex/equil-mdsim/05_npt/05_npt.gro",
        cpt=run_path+"/complex/equil-mdsim/05_npt/05_npt.cpt",
        mdp=run_path+"/complex/equil-mdsim/06_npt/06_npt.mdp",
    params:
        run_dir=run_path+"/complex/equil-mdsim/06_npt"
    output:
        gro=run_path+"/complex/equil-mdsim/06_npt/06_npt.gro",
        cpt=run_path+"/complex/equil-mdsim/06_npt/06_npt.cpt"
    threads: threads
    retries: num_retries
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

rule equil_complex_prod:
    input:
        top=input_path+"/complex/complex.top",
        ndx=input_path+"/complex/index.ndx",
        gro=run_path+"/complex/equil-mdsim/06_npt/06_npt.gro",
        cpt=run_path+"/complex/equil-mdsim/06_npt/06_npt.cpt",
        mdp=run_path+"/complex/equil-mdsim/prod/prod.mdp",
    params:
        run_dir=run_path+"/complex/equil-mdsim/prod",
    output:
        gro=run_path+"/complex/equil-mdsim/prod/prod.gro",
        cpt=run_path+"/complex/equil-mdsim/prod/prod.cpt",
        tpr=run_path+"/complex/equil-mdsim/prod/prod.tpr",
        xtc=run_path+"/complex/equil-mdsim/prod/prod.xtc",
    threads: threads
    retries: num_retries
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

rule equil_run_complex_trjconv:
    input:
        tpr=run_path+"/complex/equil-mdsim/prod/prod.tpr",
        xtc=run_path+"/complex/equil-mdsim/prod/prod.xtc"
    params:
        run_dir=run_path+"/complex/equil-mdsim/boreschcalc/"
    output:
        xtc=run_path+"/complex/equil-mdsim/boreschcalc/prod_center.xtc"
    shell:
        '''
            mkdir -p {params.run_dir}
            echo 0 | gmx trjconv -s {input.tpr} -f {input.xtc} -o {params.run_dir}/whole.xtc -pbc whole
            echo 0 | gmx trjconv -s {input.tpr} -f {params.run_dir}/whole.xtc -o {params.run_dir}/nojump.xtc -pbc nojump
            gmx trjconv -s {input.tpr} -f {params.run_dir}/nojump.xtc -o {output.xtc} -pbc mol -center -ur compact << EOF
            1
            0
            EOF
            rm {params.run_dir}/whole.xtc {params.run_dir}/nojump.xtc
        '''
        # TODO: the rm is not working

rule equil_run_complex_get_boresch_restraints:
    input:
        tpr=run_path+"/complex/equil-mdsim/prod/prod.tpr",
        xtc=run_path+"/complex/equil-mdsim/boreschcalc/prod_center.xtc"
    params:
        run_dir=run_path+"/complex/equil-mdsim/boreschcalc/",
    output:
        gro=run_path+"/complex/equil-mdsim/boreschcalc/ClosestRestraintFrame.gro",
        top=run_path+"/complex/equil-mdsim/boreschcalc/BoreschRestraint.top",
        boresch_dG_off=run_path+"/complex/equil-mdsim/boreschcalc/dG_off.dat"
    run:
        boresch.gen_restraint(
            topology = input.tpr,
            trajectory = input.xtc,
            outpath = params.run_dir,
        )
