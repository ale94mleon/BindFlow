from abfe import scripts

run_path = config["run_path"]
vdw_lambdas = config['lambdas']['ligand']['vdw']
coul_lambdas = config['lambdas']['ligand']['coul']


# Ana
rule fep_ana_gather_ligand_xvg:
    input:
        xvg_vdw_loc=expand(run_path+"/ligand/fep/simulation/vdw.{state}/prod/prod.xvg", state=range(len(config['lambdas']['ligand']['vdw']))),
        xvg_coul_loc=expand(run_path+"/ligand/fep/simulation/coul.{state}/prod/prod.xvg", state=range(len(config['lambdas']['ligand']['coul'])))
    params:
        sim_loc=run_path+"/ligand/fep/simulation",
        ana_loc=run_path+"/ligand/fep/ana",
        vdw_max_windows=len(config['lambdas']['ligand']['vdw']),
        coul_max_windows=len(config['lambdas']['ligand']['coul'])
    output:
        xvg_dir=directory(run_path+"/ligand/fep/ana/xvgs")
    shell:
        '''
            mkdir -p {params.ana_loc}/xvgs/vdw-xvg
            mkdir -p {params.ana_loc}/xvgs/coul-xvg

            # vdw
            let max_window={params.vdw_max_windows}
            for i in $(seq 0 $((max_window-1)))
            do
                cp {params.sim_loc}/vdw.${{i}}/prod/prod.xvg {params.ana_loc}/xvgs/vdw-xvg/dhdl.${{i}}.xvg
            done

            # coul
            let max_window={params.coul_max_windows}
            for i in $(seq 0 $((max_window-1)))
            do
                cp {params.sim_loc}/coul.${{i}}/prod/prod.xvg {params.ana_loc}/xvgs/coul-xvg/dhdl.${{i}}.xvg
            done
        '''

rule fep_ana_get_dg_ligand:
    input:
        xvg_dir=run_path+"/ligand/fep/ana/xvgs"
    params:
        conf_path = run_path+"/snake_conf.json",
        out_dir=run_path+"/ligand/fep/ana",
        script_dir=scripts.root_path
    output:
        complex_var=run_path+"/ligand/fep/ana/dg_results.csv"
    shell:
        "python {params.script_dir}/free_energy/calculate_ABFE_transformation_dG.py --xvg_path {input.xvg_dir} --conf_path {params.conf_path} --out_path {params.out_dir} --system_name ligand"