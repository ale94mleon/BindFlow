
input_path = config['input_data_path']
run_path = config["run_path"]
num_sim_threads = config['num_sim_threads']
num_retries = config['num_retries']
load_dependencies = config['job_extra_directives']
mdrun_extra = config['mdrun_extra_directives']

rule equil_run_complex_emin:
    input:
        top=input_path+"/complex/complex.top",
        gro=input_path+"/complex/complex.gro",
        mdp=run_path+"/complex/equil-mdsim/emin/emin.mdp"
    params:
        nthreads=num_sim_threads,
        run_dir=run_path+"/complex/equil-mdsim/emin"
    output:
        gro=run_path+"/complex/equil-mdsim/emin/emin.gro"
    threads: num_sim_threads
    retries: num_retries
    run:
        gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = input.gro,
            nthreads = params.nthreads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            cpi = True,
            **mdrun_extra
        )

rule equil_run_complex_nvt_heat:
    input:
        top=input_path+"/complex/complex.top",
        gro=run_path+"/complex/equil-mdsim/emin/emin.gro"
        mdp=run_path+"/complex/equil-mdsim/nvt_heat/nvt_heat.mdp",
    params:
        nthreads=num_sim_threads,
        run_dir=run_path+"/complex/equil-mdsim/nvt_heat"
    output:
        gro=run_path+"/complex/equil-mdsim/nvt_heat/nvt_heat.gro",
        cpt=run_path+"/complex/equil-mdsim/nvt_heat/nvt_heat.cpt"
    threads: num_sim_threads
    retries: num_retries
    run:
        gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = input.gro,
            nthreads = params.nthreads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            cpi = True,
            **mdrun_extra
        )

rule equil_run_complex_npt_eq1:
    input:
        top=input_path+"/complex/complex.top",
        gro=run_path+"/complex/equil-mdsim/nvt_heat/nvt_heat.gro",
        cpt=run_path+"/complex/equil-mdsim/nvt_heat/nvt_heat.cpt",
        mdp=run_path+"/complex/equil-mdsim/npt_equil1/npt_equil1.mdp",
    params:
        nthreads=num_sim_threads,
        run_dir=run_path+"/complex/equil-mdsim/npt_equil1"
    output:
        gro=run_path+"/complex/equil-mdsim/npt_equil1/npt_equil1.gro",
        cpt=run_path+"/complex/equil-mdsim/npt_equil1/npt_equil1.cpt"
    threads: num_sim_threads
    retries: num_retries
    run:
        gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = input.gro,
            checkpoint = input.cpt,
            nthreads = params.nthreads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            cpi = True,
            **mdrun_extra
        )

rule equil_run_complex_npt_eq2:
    input:
        top=input_path+"/complex/complex.top",
        gro=run_path+"/complex/equil-mdsim/npt_equil1/npt_equil1.gro",
        cpt=run_path+"/complex/equil-mdsim/npt_equil1/npt_equil1.cpt",
        mdp=run_path+"/complex/equil-mdsim/npt_equil2/npt_equil2.mdp",
    params:
        nthreads=num_sim_threads,
        run_dir=run_path+"/complex/equil-mdsim/npt_equil2"
    output:
        gro=run_path+"/complex/equil-mdsim/npt_equil2/npt_equil2.gro",
        cpt=run_path+"/complex/equil-mdsim/npt_equil2/npt_equil2.cpt"
    threads: num_sim_threads
    retries: num_retries
    run:
        gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = input.gro,
            checkpoint = input.cpt,
            nthreads = params.nthreads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            cpi = True,
            **mdrun_extra
        )

rule equil_run_complex_prod:
    input:
        top=run_path+"/complex/topology/complex.top",
        gro=run_path+"/complex/equil-mdsim/npt_equil2/npt_equil2.gro",
        cpt=run_path+"/complex/equil-mdsim/npt_equil2/npt_equil2.cpt"
        mdp=run_path+"/complex/equil-mdsim/npt_prod/npt_prod.mdp",
    params:
        nthreads=num_sim_threads,
        run_dir=run_path+"/complex/equil-mdsim/npt_prod",
        gmx_template=gromacs_cont_script
    output:
        gro=run_path+"/complex/equil-mdsim/npt_prod/npt_prod.gro",
        cpt=run_path+"/complex/equil-mdsim/npt_prod/npt_prod.cpt",
        tpr=run_path+"/complex/equil-mdsim/npt_prod/npt_prod.tpr",
        xtc=run_path+"/complex/equil-mdsim/npt_prod/npt_prod.xtc",
    threads: num_sim_threads
    retries: num_retries
    run:
        gmx_runner(
            mdp = input.mdp,
            topology = input.top,
            structure = input.gro,
            checkpoint = input.cpt,
            nthreads = params.nthreads,
            load_dependencies = load_dependencies,
            run_dir = params.run_dir,
            cpi = True,
            **mdrun_extra
        )

rule equil_run_complex_trjconv:
    input:
        tpr=run_path+"/complex/equil-mdsim/npt_prod/npt_prod.tpr",
        xtc=run_path+"/complex/equil-mdsim/npt_prod/npt_prod.xtc"
    params:
        run_dir=run_path+"/complex/equil-mdsim/boreschcalc/"
    output:
        xtc=run_path+"/complex/equil-mdsim/boreschcalc/npt_prod_center.xtc"
    shell:
        '''
            mkdir {params.run_dir}
            echo 0 | gmx trjconv -s {input.tpr} -f {input.xtc} -o {params.run_dir}/whole.xtc -pbc whole
            echo 0 | gmx trjconv -s {input.tpr} -f {params.run_dir}/whole.xtc -o {params.run_dir}/nojump.xtc -pbc nojump
            gmx trjconv -s {input.tpr} -f {params.run_dir}/nojump.xtc -o {output.xtc} -pbc mol -center -ur compact << EOF
            1
            0
            EOF
            rm {params.run_dir}/whole.xtc {params.run_dir}/nojump.xtc
        '''

rule equil_run_complex_get_boresch_restraints:
    input:
        tpr=run_path+"/complex/equil-mdsim/npt_prod/npt_prod.tpr",
        xtc=run_path+"/complex/equil-mdsim/boreschcalc/npt_prod_center.xtc"
    params:
        run_dir=run_path+"/complex/equil-mdsim/boreschcalc/",
        code_path = scripts.root_path
    output:
        gro=run_path+"/complex/equil-mdsim/boreschcalc/ClosestRestraintFrame.gro",
        top=run_path+"/complex/equil-mdsim/boreschcalc/BoreschRestraint.top",
        boresch_dG_off=run_path+"/complex/equil-mdsim/boreschcalc/dG_off.dat"
    shell:
        '''
            abfe-gen_boresch_restraints --top {input.tpr} --trj {input.xtc} --outpath {params.run_dir}
        '''
