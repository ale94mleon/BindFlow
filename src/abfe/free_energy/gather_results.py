import glob
import os
from abfe.utils.tools import PathLike
from typing import List
import pandas as pd
from scipy import stats
import numpy as np
from uncertainties import ufloat

def get_stats(replica_paths:List[PathLike]) -> dict:
    """Takes all the replica path and extract free energy statistics

    Parameters
    ----------
    replica_paths : List[PathLike]
        A list with all the replica paths

    Returns
    -------
    dict
        A dictionary with keywords:
        
            #. <estimator>: the value of the estimator
            #. <estimator>_sem: Standard error of the mean
            #. <estimator>_uncertainty_propagation: Propagate the uncertainties after the average.
            This use the estimated uncertainties from alchemy (Check :meth:`abfe.free_energy.analysis.run_alchemlyb`)
            #. <estimator>_num_replicas: The number of replicas employed.

    """
    # Get for each estimator the corresponded values and standard deviations
    estimator_result = dict()
    for replica_path in replica_paths:
        df = pd.read_csv(replica_path, index_col=0)
        for estimator in df.columns:
            if estimator in estimator_result:
                estimator_result[estimator].append(ufloat(df.loc['value', estimator],df.loc['std_dev', estimator]))
            else:
                estimator_result[estimator] = [ufloat(df.loc['value', estimator],df.loc['std_dev', estimator])]
    
    # Build the final result
    final_result = dict()
    for estimator in estimator_result:
        mean = np.mean(estimator_result[estimator])
        # Save results
        final_result[estimator] = mean.nominal_value
        final_result[f"{estimator}_sem"] = stats.sem([value.nominal_value for value in estimator_result[estimator]], ddof=1)
        final_result[f"{estimator}_uncertainty_propagation"] = mean.std_dev
        final_result[f"{estimator}_num_replicas"] = len(estimator_result[estimator])
    
    return final_result

def get_all_dgs(root_folder_path:PathLike, out_csv:PathLike = None) -> pd.DataFrame:
    """Get the independent free energy results and gather them

    Parameters
    ----------
    root_folder_path : PathLike
        Where the simulation run. Inside of it should be the files: root_folder_path + "/*/*/dG_results.csv".
        THis directory is the same specified on :meth:`abfe.calculate_abfe.calculate_abfe` through the keyword `out_root_folder_path`
    out_csv : PathLike, optional
        If given a pandas.DataFrame will be written as csv file, by default None

    Returns
    -------
    pd.DataFrame
        All gather results. In case there are not any dG_results.csv. It will return None
    """
    # Get all dG_results.csv files
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
            statistics = get_stats(dg_files_dir[ligand])
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
if __name__ == '__main__':
    pass
    root_folder_path = "/home/uds_alma015/GIT/BindFlow/examples/internal_example/CyclophilinD/abfe"
    get_all_dgs(root_folder_path)