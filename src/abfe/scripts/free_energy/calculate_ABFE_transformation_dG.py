"""
calculate dF for complex system
"""

import argparse
import json

import numpy as np

from abfe.scripts.free_energy.alchemlyb_analysis import analyze_ligand


def calculate_transformation_dG(conf_path: str, xvg_path: str,
                                system_name: str,
                                boresch_data: str, out_path: str,
                                min_samples: int = 200, lower_bound_ndatapoints: int = 1000, ):
    conf = json.load(open(conf_path, "r"))

    if (system_name == "complex"):
        system_steps_windows = {
            "bonded-xvg": conf['n_rest_windows_complex'],
            'vdw-xvg': conf['n_vdw_windows_complex'],
            'coul-xvg': conf['n_coul_windows_complex']
        }
    elif (system_name == "ligand"):
        system_steps_windows = {
            'vdw-xvg': conf['n_vdw_windows_ligand'],
            'coul-xvg': conf['n_coul_windows_ligand']
        }

    else:
        raise ValueError("The provided system_name is unknown, please provide complex or ligand")

    res_df = analyze_ligand(prefix=xvg_path,
                            system_steps_windows=system_steps_windows,
                            system_name=system_name,
                            lower=lower_bound_ndatapoints,
                            min_samples=min_samples)

    # include boresch correction
    print(boresch_data)
    if (boresch_data is not None):
        col = []
        res_V = float(np.loadtxt(boresch_data))
        for row in res_df.index:
            if (row == "sys"):
                col.append("ligand")
            elif (row == "windows"):
                col.append("-")
            else:
                col.append(res_V)
        res_df['boresch'] = col

    res_df.to_csv(out_path + "/dg_results.csv")
    return out_path


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--xvg_path', default='./',
                        help='input xvg path')
    parser.add_argument('--conf_path', default='../../../snake_conf.json',
                        help='path to job conf path')
    parser.add_argument('--out_path', default='./',
                        help='output path for writing files')
    parser.add_argument('--system_name', default='test',
                        help='output path for writing files')
    parser.add_argument("--boresch_data",
                        help='boresch restraint correction', required=False, default=None)
    parser.add_argument("--lower_bound_ndatapoints",
                        help='minimal traj length', required=False, default=1000)
    parser.add_argument("--min_samples",
                        help='minimal number samples', required=False, default=200)

    args = parser.parse_args()

    calculate_transformation_dG(conf_path=args.conf_path,
                                xvg_path=args.xvg_path,
                                out_path=args.out_path,
                                boresch_data=args.boresch_data,
                                system_name=args.system_name,
                                min_samples=args.min_samples,
                                lower_bound_ndatapoints=args.lower_bound_ndatapoints
                                )
