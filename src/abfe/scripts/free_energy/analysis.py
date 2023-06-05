"""
calculate dF for ligand system
"""
import math
import os
import warnings
import json
import numpy as np
import pandas as pd
from uncertainties import ufloat
from alchemlyb import concat
from alchemlyb.estimators import TI, MBAR
from alchemlyb.parsing.gmx import extract_dHdl, extract_u_nk
from alchemlyb.preprocessing import statistical_inefficiency, slicing
from alchemlyb.visualisation import plot_convergence
from alchemlyb.convergence import forward_backward_convergence
from alchemlyb.postprocessors import units
from abfe.utils.tools import PathLike

def run_alchemlyb(xvgs: list, lower: int = None, upper: int = None, min_samples: int = 500, temperature: float = 298.15, convergency_plots_prefix:str = None):
    """
    Function to get MBAR and TI estimates using alchemlyb from an input set of
    xvgs

    Parameters
    ----------
    xvgs : list of str
        list of filenames for input xvg files. This list must be sorted based on the lambda values.
        For example. If we have coul and vdw and three lambda points; the list must be:
        either: [coul1, coul2, coul3, vdw1, vdw2, vdw3] or [vdw1, vdw2, vdw3, coul1, coul2, coul3]
        TODO: We can sort the input files too. I am also not sure that that is needed. Maybe alchemlyb 
        does the sort internally.
    overlap_path : str
        path to write overlap matrix (if None, no matrix will be written)
        [None]
    lower : int
        starting time to sample dhdl xvgs from [None]
    upper : int
        inclusive end time to sample dhdl xvgs from [None]
    min_samples : int
        minimum number of samples to analyze, if statistical inefficiency
        returns fewer samples than this, then min_samples will be picked
        instead [500]
    temperature : float
        simulation temperature [298.15]
    convergency_plots_prefix : str
        If gives, it will plot the convergency to {convergency_plots_prefix}convergence_TI.pdf and {convergency_plots_prefix}convergence_MBAR.pdf, by default None

    Returns
    -------
    deltaG : dict
        deltaG = {
            'MBAR': {
                'value': <value>,
                'error': <error>
            },
            'TI': {
                'value': <value>,
                'error': <error>
            },
        }
        two entry dictionary containing the MBAR and TI free energy and
        associated variance error estimate in kcal/mol
    """

    # Check how many independent samples do we have in our data.
    sub_steps = []
    for xvg in xvgs:
        extracted_dHdls = extract_dHdl(xvg, T = temperature)
        
        df = slicing(extracted_dHdls, lower = lower, upper = upper)
        
        # Calculate the time step on the data based on statistical_inefficiency (or autocorrelation)
        df_ineff = statistical_inefficiency(df, series=df.iloc[:, 0])
        if len(df_ineff) != 0:
            ineff_step = math.ceil(len(df) / len(df_ineff))
        else:
            ineff_step = 1
            wmsg = "statistical_inefficiency does not give data. This usually means that the sampling is too poor.\n"\
                          f"xvg: {xvg}"
            warnings.warn(wmsg)

        
        # Check the lag of the step to fulfill the minimum number of samples
        step_cutoff = int(len(df) / min_samples)
        if step_cutoff == 0:
            step_cutoff = 1
            wmsg = f"The number of raw data point ({len(df)}) is less than min_samples = {min_samples}. This usually means that the sampling is too poor\n"\
                          f"xvg: {xvg}"
            warnings.warn(wmsg)

        # Select the proper step lag
        if ineff_step <= step_cutoff:
            sub_steps.append(ineff_step)
        else:
            sub_steps.append(int(step_cutoff))

    print(f"number of samples per window: {[int(len(df) / sub_step) for sub_step in sub_steps]}")

    dhdls_data = [slicing(extract_dHdl(xvg, T=temperature), lower=lower, upper=upper, step=step) for xvg, step in zip(xvgs, sub_steps)]
    u_nks_data = [slicing(extract_u_nk(xvg, T=temperature), lower=lower, upper=upper, step=step) for xvg, step in zip(xvgs, sub_steps)]

    # Get estimations
    mbar = MBAR(maximum_iterations=1000000).fit(concat(u_nks_data))
    ti = TI().fit(concat(dhdls_data))

    # TODO And print some images that show some convergence related factors
    # Convert values and errors to kcal/mol, and access the free energy difference between the states at lambda 0.0 and 1.0
    deltaG = {
        'MBAR': {
            'value': units.to_kcalmol(mbar.delta_f_).iloc[0, -1],
            'error': units.to_kcalmol(mbar.d_delta_f_).iloc[0, -1]
        },
        'TI': {
            'value': units.to_kcalmol(ti.delta_f_).iloc[0, -1],
            'error': units.to_kcalmol(ti.d_delta_f_).iloc[0, -1]
        }       
    }
    # # Evaluate convergency
    # # TODO: On MBAR I am getting LLVM ERROR: pthread_create failed: Resource temporarily unavailable Aborted
    # # And I can not capture this error
    # if convergency_plots_prefix:
    #     try:
    #         ax = plot_convergence(forward_backward_convergence(dhdls_data, 'TI'))
    #         ax.figure.savefig(f'{convergency_plots_prefix}convergence_TI.pdf')
    #         # ax = plot_convergence(forward_backward_convergence(u_nks_data, 'MBAR'))
    #         # ax.figure.savefig(f'{convergency_plots_prefix}convergence_MBAR.pdf')
    #     except Exception as e:
    #         print(f"Not possible to evaluate convergency on: {xvgs}\n. Exception {e} was got it")
    return deltaG

# Used on complex(ligand)_fep_ana.smk
def get_dG_contributions(
        boresch_data:PathLike = None,
        out_json_path:PathLike = 'dg_contributions.json',
        lower: int = None,
        upper: int = None,
        min_samples: int = 500,
        temperature: float = 298.15,
        convergency_plots_prefix:str = None,
        **kwargs):
    """It calculate and gather the vdw, coul and bonded (in case of a complex) dG' contributions.
    It also adds the analytical correction due to the ligand restraints.

    Parameters
    ----------
    boresch_data : PathLike, optional
        File with the boresch analytical corrections (this should only be used for the ligand 
        for the ligand), by default None
    out_csv_path : PathLike, optional
        Path to output the results, by default 'dg_contributions.csv'. It is going to have as columns: [MBAR, TI, boresch]
        and as index: [value, error].
    lower : int, optional
         Upper time to slice, by default None
    upper : int, optional
        Upper time to slice to (inclusive), by default None
    min_samples : int, optional
        Minimum number of samples to use, by default 500
    temperature : float, optional
        Temperature of the simulation, by default 298.15
    convergency_plots_prefix : str
        If gives, it will plot the convergency to {convergency_plots_prefix}{<kwargs_name>}_convergence_TI.pdf and {convergency_plots_prefix}{<kwargs_name>}_convergence_MBAR.pdf, by default None
    **kwargs : optional
        Only vdw = <List[PathLike]>, coul = <List[PathLike]> and bonded = <List[PathLike]> are valid extra keywords.
        Those are the path to the xvg files of the corresponded lambda type.
        # TODO: They should be order form lambda 0 to 1??. it is safe to do it.
    Raises
    ------
    ValueError
        Invalid kwargs. Only valid: vdw, could and bonded
    ValueError
        The value of the kwargs is not a list
    FileNotFoundError
        If some xvg files are not found.
    """
    # Check validity of keywords:
    valid_kwargs = ['vdw', 'coul','bonded']
    for key in kwargs:
        if key not in valid_kwargs:
            raise ValueError(f"The provided extra keyword '{key}' is not valid. Choose from: {valid_kwargs}")
        elif not isinstance(kwargs[key], list):
            raise ValueError(f"The provided extra keyword '{key}' is valid. But its value '{kwargs[key]}' is not an integer")

    system_results = dict()
    # For each lambda_type get each its xvg in a sorted manner.
    for lambda_type in kwargs:
        for xvg_file in kwargs[lambda_type]:
            # Check that all xvg files exist
            if not os.path.isfile(xvg_file):
                raise FileNotFoundError(f"Provided xvg file: {xvg_file} for lambda_type = {lambda_type} ")
        if convergency_plots_prefix:
            convergency_plots_prefix_to_use = f"{convergency_plots_prefix}{lambda_type}_"
        dG = run_alchemlyb(xvgs = kwargs[lambda_type], lower=lower, upper=upper, min_samples=min_samples, temperature=temperature, convergency_plots_prefix=convergency_plots_prefix_to_use)
        
        ddG_estimator = abs(dG['MBAR']['value'] - dG['TI']['value'])
        if ddG_estimator > 0.5:
            wmsg = (f'|dG_MBAR - dG_TI| = {ddG_estimator} > 0.5 kcal/mol for {lambda_type}')
            warnings.warn(wmsg)

        system_results[lambda_type] = dG

    # include boresch correction
    if boresch_data:
        system_results['boresch'] = float(np.loadtxt(boresch_data))
    # Write the data
    with open(out_json_path, 'w') as out:
        json.dump(system_results, out, indent=4)


def get_ufloat(value_error:dict) -> ufloat:
    """Simple function to get the ufloat class
    form a dict with keywords value and error

    Parameters
    ----------
    value_error : dict
        The dict with keywords: value and error

    Returns
    -------
    ufloat
        The ufloat class ready to propagate uncertainties.
    """
    return ufloat(value_error['value'], value_error['error'])

# TODO, check what is going on here calculate_ABFE_ligand_dG
def get_dg_cycle(ligand_contributions:PathLike = 'dg_ligand_contributions.json', complex_contributions:PathLike = 'dg_complex_contributions.json', out_csv:PathLike = 'dG_results.csv'):

    # Create new dict containing the results and calculate complete process:
    with open(ligand_contributions, 'r') as lc:
        ligand_dict = json.load(lc)
    with open(complex_contributions, 'r') as cc:
        complex_dict = json.load(cc)

    # TODO WARNING!! I am not completely sure about the sign of boresch, and to be honest more or less for the rest either
    # TODO, I can do this in a more dynamic way
    dG_MBAR = get_ufloat(ligand_dict["coul"]["MBAR"]) + get_ufloat(ligand_dict["vdw"]["MBAR"]) + ligand_dict["boresch"] + \
        get_ufloat(complex_dict["vdw"]["MBAR"]) + get_ufloat(complex_dict["coul"]["MBAR"]) + get_ufloat(complex_dict["bonded"]["MBAR"])

    dG_TI = get_ufloat(ligand_dict["coul"]["TI"]) + get_ufloat(ligand_dict["vdw"]["TI"]) + ligand_dict["boresch"] + \
        get_ufloat(complex_dict["vdw"]["TI"]) + get_ufloat(complex_dict["coul"]["TI"]) + get_ufloat(complex_dict["bonded"]["TI"])
    
    deltaG = {
        'MBAR': {
            'value': dG_MBAR.nominal_value,
            'std_dev': dG_MBAR.std_dev
        },
        'TI': {
            'value': dG_TI.nominal_value,
            'std_dev': dG_TI.std_dev
        }       
    }
    pd.DataFrame(deltaG).to_csv(out_csv)

if __name__ == "__main__":
    pass
    # import glob
    # from abfe.mdp import mdp
    
    # # Parameters complex
    # run_path = "/scratch/uds_alma015/smaug/data/users/alejandro/simulation/BindFlow_simulations/biotin-streptavidin/abfe/"
    # mdp_vdw_0_prod = run_path + "/biotin_a/1/complex/fep/simulation/vdw.0/prod/prod.mdp"

    # xvg_vdw_loc = [os.path.join(path, "prod/prod.xvg") for path in glob.glob(run_path+"/biotin_a/1/complex/fep/simulation/vdw.*")]
    # coul_vdw_loc = [os.path.join(path, "prod/prod.xvg") for path in glob.glob(run_path+"/biotin_a/1/complex/fep/simulation/coul.*")]
    # bonded_vdw_loc = [os.path.join(path, "prod/prod.xvg") for path in glob.glob(run_path+"/biotin_a/1/complex/fep/simulation/bonded.*")]


    # boresch_data = None
    # # Get the simulaiton temperature from the prod.mdp of the state 0 of vdw
    # temperature = float(mdp.MDP().from_file(mdp_vdw_0_prod).parameters['ref-t'].split()[0]) 
    # get_dG_contributions(

    #     boresch_data = boresch_data,
    #     out_json_path = 'results/dg_complex_contributions.json',
    #     # Check if it is necessary to remove some initial burning simulation time
    #     lower = 4000,
    #     upper = None,
    #     min_samples = 500,
    #     temperature = temperature,
    #     convergency_plots_prefix = 'results/complex_',
    #     vdw = sorted(xvg_vdw_loc, key=lambda x: int(os.path.normpath(x).split(os.path.sep)[-3].split('.')[-1])),
    #     coul = sorted(coul_vdw_loc, key=lambda x: int(os.path.normpath(x).split(os.path.sep)[-3].split('.')[-1])),
    #     bonded = sorted(bonded_vdw_loc, key=lambda x: int(os.path.normpath(x).split(os.path.sep)[-3].split('.')[-1])),
    # )


    # # Parameters ligand
    # mdp_vdw_0_prod = run_path + "/biotin_a/1/ligand/fep/simulation/vdw.0/prod/prod.mdp"

    # xvg_vdw_loc = [os.path.join(path, "prod/prod.xvg") for path in glob.glob(run_path+"/biotin_a/1/ligand/fep/simulation/vdw.*")]
    # coul_vdw_loc = [os.path.join(path, "prod/prod.xvg") for path in glob.glob(run_path+"/biotin_a/1/ligand/fep/simulation/coul.*")]

    # boresch_data = run_path + "/biotin_a/1/complex/equil-mdsim/boreschcalc/dG_off.dat"
    # # Get the simulaiton temperature from the prod.mdp of the state 0 of vdw
    # temperature = float(mdp.MDP().from_file(mdp_vdw_0_prod).parameters['ref-t'].split()[0]) 
    # get_dG_contributions(
    #     boresch_data = boresch_data,
    #     out_json_path = 'results/dg_ligand_contributions.json',
    #     # Check if it is necessary to remove some initial burning simulation time
    #     lower = 4000,
    #     upper = None,
    #     min_samples = 500,
    #     temperature = temperature,
    #     convergency_plots_prefix = 'results/ligand_',
    #     vdw = sorted(xvg_vdw_loc, key=lambda x: int(os.path.normpath(x).split(os.path.sep)[-3].split('.')[-1])),
    #     coul = sorted(coul_vdw_loc, key=lambda x: int(os.path.normpath(x).split(os.path.sep)[-3].split('.')[-1])),
    # )


    # get_dg_cycle(
    #     ligand_contributions='results/dg_ligand_contributions.json',
    #     complex_contributions='results/dg_complex_contributions.json',
    #     out_csv='results/dG_results.csv'
    # )