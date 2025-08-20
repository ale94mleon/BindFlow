import glob
import json
import os
from pathlib import Path
from typing import List, Union

import numpy as np
import pandas as pd

from bindflow.utils.tools import PathLike, sum_uncertainty_propagation


def get_fep_stats(replica_paths: List[PathLike]) -> dict:
    """Takes all the replica path and extract free energy statistics

    Parameters
    ----------
    replica_paths : List[PathLike]
        A list with all the replica paths

    Returns
    -------
    dict
        A dictionary with keywords:

            #. <estimator>: the average value of the estimator
            #. <estimator>_sem: Standard error of the mean
            #. <estimator>_uncertainty_propagation: Propagate the uncertainties after the average.
            This use the estimated uncertainties from alchemlyb (Check :func:`bindflow.free_energy.fep_analysis.run_alchemlyb`)
            #. <estimator>_num_replicas: The number of replicas employed.

    """
    # Get for each estimator the corresponded values and standard deviations
    estimator_result = dict()
    for replica_path in replica_paths:
        df = pd.read_csv(replica_path, index_col=0)
        for estimator in df.columns:
            if estimator in estimator_result:
                estimator_result[estimator].append((df.loc['value', estimator], df.loc['std_dev', estimator]))
            else:
                estimator_result[estimator] = [(df.loc['value', estimator], df.loc['std_dev', estimator])]

    # Build the final result
    final_result = dict()
    for estimator in estimator_result:
        mean_value = np.mean([value_error[0] for value_error in estimator_result[estimator]])
        mean_std_dev = sum_uncertainty_propagation(
            errors=[value_error[1] for value_error in estimator_result[estimator]],
            coefficients=len(estimator_result[estimator]) * [1/len(estimator_result[estimator])]
        )
        # Save results
        final_result[estimator] = mean_value
        final_result[f"{estimator}_sem"] = pd.Series([value_error[0] for value_error in estimator_result[estimator]]).sem(ddof=1)
        final_result[f"{estimator}_uncertainty_propagation"] = mean_std_dev
        final_result[f"{estimator}_num_replicas"] = len(estimator_result[estimator])

    return final_result


def get_all_fep_dgs(root_folder_path: PathLike, out_csv: PathLike = None) -> pd.DataFrame:
    """Get the independent FEP results and gather them.
    Average and standard error of the mean are reported.

    Parameters
    ----------
    root_folder_path : PathLike
        Where the simulation run. Inside it should be the files: root_folder_path + "/*/*/dG_results.csv".
        This directory is the same specified on :func:`bindflow.runnner.calculate` through the keyword `out_root_folder_path`
    out_csv : PathLike, optional
        If given a pandas.DataFrame will be written as csv file, by default None

    Returns
    -------
    pd.DataFrame
        All gather results. If there are not dG_results.csv; It will return an empty DataFrame
    """
    # Get all dG_results.csv files
    root_folder_path = str(root_folder_path)
    dg_files_dir = dict()
    for dg_file in glob.glob(root_folder_path + "/*/*/dG_results.csv"):
        dg_file = os.path.normpath(dg_file)
        ligand_name = dg_file.split(os.path.sep)[-3]
        if ligand_name in dg_files_dir:
            dg_files_dir[ligand_name].append(dg_file)
        else:
            dg_files_dir[ligand_name] = [dg_file]

    gathered_results = []
    if dg_files_dir:
        for ligand in dg_files_dir:
            statistics = get_fep_stats(dg_files_dir[ligand])
            statistics['ligand'] = ligand
            gathered_results.append(statistics)

        gathered_results = pd.DataFrame(gathered_results)
        # Put the column 'ligand' at the beginning
        columns = ['ligand'] + [col for col in gathered_results.columns if col != 'ligand']
        gathered_results = gathered_results[columns]
        # Safe data on request
        if out_csv:
            gathered_results.to_csv(out_csv)
        return gathered_results
    else:
        print(f"There is not dG_results.csv yet on {root_folder_path}/*/*")
        return pd.DataFrame()


def get_raw_fep_data(root_folder_path: PathLike, out_csv: PathLike = None) -> pd.DataFrame:
    """Generate raw dat for an FEP calculation using BindFlow

    Parameters
    ----------
    root_folder_path : PathLike
        Where the simulation run. Inside it should be the files: root_folder_path + "/*/*/complex/fep/ana/dg_complex_contributions.json".
        This directory is the same specified on :func:`bindflow.runners.calculate` through the keyword `out_root_folder_path`
    out_csv : PathLike, optional
        If given a pandas.DataFrame will be written as csv file, by default None

    Returns
    -------
     pd.DataFrame
        Raw FEP data, all the contributions for all ligand/replicas
    """
    sample_data = []
    root_folder_path = Path(root_folder_path).resolve()
    for item1 in root_folder_path.iterdir():
        if item1.is_dir():
            ligand = item1.stem
            for item2 in item1.iterdir():
                if item2.is_dir():
                    replica = item2.stem
                    complex_json = item2/"complex/fep/ana/dg_complex_contributions.json"
                    ligand_json = item2/"ligand/fep/ana/dg_ligand_contributions.json"
                    if complex_json.is_file() and ligand_json.is_file():
                        with open(complex_json, 'r') as cj:
                            complex_data = json.load(cj)
                        with open(ligand_json, 'r') as lj:
                            ligand_data = json.load(lj)

                        sample_data.append(
                            [
                                ligand,
                                replica,
                                complex_data['vdw']['MBAR']['value'],
                                complex_data['coul']['MBAR']['value'],
                                complex_data['bonded']['MBAR']['value'],
                                ligand_data['vdw']['MBAR']['value'],
                                ligand_data['coul']['MBAR']['value'],

                                complex_data['vdw']['TI']['value'],
                                complex_data['coul']['TI']['value'],
                                complex_data['bonded']['TI']['value'],
                                ligand_data['vdw']['TI']['value'],
                                ligand_data['coul']['TI']['value'],

                                complex_data['boresch'],

                                complex_data['vdw']['MBAR']['error'],
                                complex_data['coul']['MBAR']['error'],
                                complex_data['bonded']['MBAR']['error'],
                                ligand_data['vdw']['MBAR']['error'],
                                ligand_data['coul']['MBAR']['error'],

                                complex_data['vdw']['TI']['error'],
                                complex_data['coul']['TI']['error'],
                                complex_data['bonded']['TI']['error'],
                                ligand_data['vdw']['TI']['error'],
                                ligand_data['coul']['TI']['error'],
                            ]
                        )
    df = pd.DataFrame(
        sample_data,
        columns=[
            'ligand', 'replica',
            'mbar_complex_vdw_value', 'mbar_complex_coul_value', 'mbar_complex_bonded_value', 'mbar_ligand_vdw_value', 'mbar_ligand_coul_value',
            'ti_complex_vdw_value', 'ti_complex_coul_value', 'ti_complex_bonded_value', 'ti_ligand_vdw_value', 'ti_ligand_coul_value',
            'boresch',
            'mbar_complex_vdw_error', 'mbar_complex_coul_error', 'mbar_complex_bonded_error', 'mbar_ligand_vdw_error', 'mbar_ligand_coul_error',
            'ti_complex_vdw_error', 'ti_complex_coul_error', 'ti_complex_bonded_error', 'ti_ligand_vdw_error', 'ti_ligand_coul_error',
        ]
        )
    if out_csv:
        df.to_csv(out_csv)
    return df


def get_all_mmxbsa_dgs(full_df: pd.DataFrame, columns_to_process: Union[None, List[str]] = None, out_csv: PathLike = None) -> pd.DataFrame:
    """Get the independent MM(P/G)BSA free energy results and gather them.
    Average and standard error of the mean are across all replicas and samples for each ligand are reported

    Parameters
    ----------
    full_df : pd.DataFrame,
        DataFrame generated by :func:`bindflow.free_energy.gather_results.get_raw_mmxbsa_dgs`
    columns_to_process :  Union[None, List[str]], optional
        The columns of full_df to process, by default None which means that the following will be used:
            "dg_c2_pb", "dg_c2_gb", "dg_ie_pb", "dg_ie_gb",
            "dg_qh_pb", "dg_qh_gb", "dg_en_pb", "dg_en_gb",
            "c2_pb", "c2_gb", "ie_pb", "ie_gb", "qh"
    out_csv : PathLike, optional
        If given a pandas.DataFrame will be written as csv file, by default None

    Returns
    -------
    pd.DataFrame
        All gather results. In case there are not dG_results.csv. It will return an empty DataFrame
    """
    if len(full_df):
        if columns_to_process is None:
            columns_to_process = [
                "dg_c2_pb",
                "dg_c2_gb",
                "dg_ie_pb",
                "dg_ie_gb",
                "dg_qh_pb",
                "dg_qh_gb",
                "dg_en_pb",
                "dg_en_gb",
                "c2_pb",
                "c2_gb",
                "ie_pb",
                "ie_gb",
                "qh"
                ]

        # # Convert replica and sample column to integers
        # full_df['sample'] = full_df['sample'].astype(int)

        # Group by 'name' and calculate mean and SEM
        grouped = full_df.groupby('name')

        mean_df = grouped[columns_to_process].mean().reset_index()
        sem_df = grouped[columns_to_process].sem().reset_index()

        # Count unique replicas
        replica_counts = grouped['replica'].nunique().reset_index(name='num_replicas')

        # Calculate total samples
        total_samples = grouped['sample'].count().reset_index(name='total_samples')

        # Merge the results
        final_df = mean_df.merge(sem_df, on='name', suffixes=('_mean', '_sem'))
        final_df = final_df.merge(replica_counts, on='name')
        final_df = final_df.merge(total_samples, on='name')

        if out_csv:
            final_df.to_csv(out_csv, index=False)
        return final_df
    else:
        return pd.DataFrame()


def get_raw_mmxbsa_dgs(root_folder_path: PathLike, out_csv: PathLike = None) -> pd.DataFrame:
    """Main function to retrieve MM(P/G)BSA simulation data
    generated through BindFlow

    Parameters
    ----------
    root_folder_path : PathLike
        Where the simulation run. Inside it should be the files: root_folder_path + "/*/*/complex/mmpbsa/simulation/*/mmxbsa.csv".
        This directory is the same specified on :func:`bindflow.run_mmpbsa.calculate_mmpbsa` through the keyword `out_root_folder_path`
    out_csv : PathLike, optional
        If given, a pandas.DataFrame will be written as csv file with the raw data, by default None

    Returns
    -------
    pd.DataFrame
        Raw MM(P/G)BSA data
    """
    from bindflow.free_energy import mmxbsa_analysis
    root_folder_path = str(root_folder_path)
    collected_dfs = []
    collected_files = glob.glob(root_folder_path + "/*/*/complex/mmpbsa/simulation/*/mmxbsa.csv")
    if len(collected_files) == 0:
        print(f"There is not mmxbsa.csv yet on {root_folder_path}/*/*/complex/mmpbsa/simulation/*/mmxbsa.csv")
        return pd.DataFrame()
    for inp_file in collected_files:  # collecting all of the mmxbsa.csv files
        string_data = inp_file.removeprefix(root_folder_path).split("/")
        ligand_name, replica, sample = string_data[1], string_data[2], string_data[6].removeprefix("rep.")
        collected_dfs.append(mmxbsa_analysis.convert_format_flatten(pd.read_csv(inp_file), ligand_name, replica, sample))
    full_df = pd.concat(collected_dfs, ignore_index=True)

    if out_csv:
        full_df.to_csv(out_csv, index=False)
    return full_df


if __name__ == '__main__':
    pass
