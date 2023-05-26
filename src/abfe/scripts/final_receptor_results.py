import argparse
import glob
import os

import pandas as pd

# TODO I am here, change this code based on the modifications that I made.
# Load all individual results:
def get_all_subresults(in_root_dir: str):
    files = glob.glob(in_root_dir + "/*/*/dG_results.csv")

    dfs = []
    for f in files:
        print(f)
        num_replicate = f.split("/")[-2]
        name = f.split("/")[-3].replace("ligand-", "")

        df = pd.read_csv(f, index_col=0)

        df['ligand'] = name
        df['replicate'] = num_replicate
        dfs.append(df)
    return dfs


def extract_final_results(df_app: pd.DataFrame):
    df_abs = df_app.where(df_app.step == "ABFE").dropna()
    ds = []
    for lig in df_abs.ligand.unique():
        tmp_df = df_abs.where(df_abs.ligand == lig).dropna()
        d = {"ligand": lig,
             "ABFE_mean": tmp_df.MBAR.mean().round(2),
             "ABFE_err": tmp_df.MBAR.std().round(2),
             "nreplicates": tmp_df.shape[0],
             }
        ds.append(d)

    df_final = pd.DataFrame(ds)
    return df_final


def get_final_results(in_root_dir: str, out_dir: str):
    if (not os.path.exists(in_root_dir)):
        raise IOError("Could not find the input directory: < " + str(in_root_dir) + " >")

    dfs = get_all_subresults(in_root_dir=in_root_dir)

    if (len(dfs) == 0):
        raise ValueError("no results were found in directory: " + str(in_root_dir))
    df_app = pd.concat(dfs, ignore_index=True)
    df_final = extract_final_results(df_app=df_app)

    out_df_final_results = out_dir + "/abfe_results.csv"
    out_df_single_detailed_results = out_dir + "/abfe_single_detailed_results.csv"

    df_final.to_csv(out_df_final_results, sep="\t", na_rep="-")
    df_app.to_csv(out_df_single_detailed_results, sep="\t", na_rep="-")

    return out_df_final_results, out_df_single_detailed_results


if (__name__ == "__main__"):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in_root_dir', required=True, help='where are the ligand folders containing abfe results')
    parser.add_argument('-o', '--out_dir', required=False, default='.',
                        help='where to output the result file')
    args = parser.parse_args()

    out_dir = os.path.abspath(args.out_dir)
    in_root_dir = os.path.abspath(args.in_root_dir)

    out_df_final_results, out_df_single_detailed_results = get_final_results(out_dir=out_dir, in_root_dir=in_root_dir)
    print("writing out: ", out_df_final_results, "\n\t", out_df_single_detailed_results)
