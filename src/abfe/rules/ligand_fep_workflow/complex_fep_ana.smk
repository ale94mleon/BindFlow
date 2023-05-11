from abfe import scripts

#Ana
rule fep_ana_gather_complex_xvg:
    input:
        xvg_vdw_loc=expand(run_path+"/complex/fep/simulation/vdw.{state}/prod/prod.xvg", state=range(len(config['lambdas']['complex']['vdw']))),
        xvg_coul_loc=expand(run_path+"/complex/fep/simulation/coul.{state}/prod/prod.xvg", state=range(len(config['lambdas']['complex']['coul']))),
        xvg_bonded_loc=expand(run_path+"/complex/fep/simulation/coul.{state}/prod/prod.xvg", state=range(len(config['lambdas']['complex']['bonded'])))

    params:
        sim_loc=run_path+"/complex/fep/simulation",
        ana_loc=run_path+"/complex/fep/ana",
        vdw_max_windows=len(config['lambdas']['complex']['vdw']),
        coul_max_windows=len(config['lambdas']['complex']['coul']),
        bonded_max_windows=len(config['lambdas']['complex']['bonded']),
    output:
        xvg_dir=directory(run_path+"/complex/fep/ana/xvgs")
    shell:
        '''
            mkdir -p {params.ana_loc}/xvgs/vdw-xvg
            mkdir -p {params.ana_loc}/xvgs/coul-xvg
            mkdir -p {params.ana_loc}/xvgs/bonded-xvg

            # vdw
            let max_window={params.vdw_max_windows}  
            for i in $(seq 0 1 $((max_window-1)))
            do
                cp {params.sim_loc}/vdw.${{i}}/prod/prod.xvg {params.ana_loc}/xvgs/vdw-xvg/dhdl.${{i}}.xvg
            done

            # coul
            let max_window={params.coul_max_windows}  
            for i in $(seq 0 1 $((max_window-1)))
            do
                cp {params.sim_loc}/coul.${{i}}/prod/prod.xvg {params.ana_loc}/xvgs/coul-xvg/dhdl.${{i}}.xvg
            done

            # bonded
            let max_window={params.bonded_max_windows}  
            for i in $(seq 0 1 $((max_window-1)))
            do
                cp {params.sim_loc}/bonded.${{i}}/prod/prod.xvg {params.ana_loc}/xvgs/bonded-xvg/dhdl.${{i}}.xvg
            done
        '''

rule fep_ana_get_dg_complex:
    input:
        xvg_dir=run_path+"/complex/fep/ana/xvgs",
        boresch_dat = run_path+"/complex/equil-mdsim/boreschcalc/dG_off.dat"
    params:
        conf_path = run_path+"/snake_conf.json",
        out_dir=run_path+"/complex/fep/ana",
        script_dir = scripts.root_path
    output:
        complex_var=run_path+"/complex/fep/ana/dg_results.csv"
    shell:
        "python {params.script_dir}/free_energy/calculate_ABFE_transformation_dG.py --xvg_path {input.xvg_dir}  --conf_path {params.conf_path} --out_path {params.out_dir} --boresch_data {input.boresch_dat} --system_name complex"

