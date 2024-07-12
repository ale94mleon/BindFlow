#!/usr/bin/env python3

import argparse
import logging
import os
from bindflow import __version__

loggers = [logging.getLogger(name) for name in logging.root.manager.loggerDict]
for logger in loggers:
    logger.setLevel(logging.NOTSET)


def dag_maker(input_path, out_name):
    from bindflow.utils import tools
    cwd = os.getcwd()
    os.chdir(input_path)
    tools.run(f"snakemake --dag | dot -Tpng -o {out_name}.png", interactive=True)
    os.chdir(cwd)


def abfe_check_results(out_root_folder_path, out_csv_summary, out_csv_raw):
    from bindflow.free_energy import gather_results

    df_summary = gather_results.get_all_abfe_dgs(root_folder_path=out_root_folder_path)
    if len(df_summary):
        df_summary = df_summary.sort_values(by='MBAR').reset_index()
        if out_csv_summary:
            df_summary.to_csv(out_csv_summary)
        print(df_summary)
    if out_csv_raw:
        df_raw = gather_results.get_raw_abfe_data(root_folder_path=out_root_folder_path)
        if len(df_raw):
            df_raw.to_csv(out_csv_raw)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--version',
        action='version',
        version=f"BindFlow: {__version__}")

    subparsers = parser.add_subparsers(required=False, dest="command")

    dag = subparsers.add_parser('dag', help="Build the DAG of the workflow")
    dag.add_argument(
        '-i',
        dest='input_path',
        help='Where should the dag command should be executed. The path where the main Snakemake file is located',
        type=str, default='.')
    dag.add_argument(
        '-o',
        dest='out_name',
        help='Name of the output image. The suffix `.png` will be added at the end',
        default='dag', type=str)

    abfe_check = subparsers.add_parser(
        'check_abfe',
        help="Check for completion of an ABFE workflow")
    abfe_check.add_argument(
        dest='out_root_folder_path',
        help='ABFE directory (`out_root_folder_path` kwarg of :meth:`bindflow.run_abfe.calculate_abfe`)',
        type=str)
    abfe_check.add_argument(
        '-os', '--out_csv_summary',
        help="The path to output the summary csv file, by default None",
        dest='out_csv_summary',
        nargs=argparse.OPTIONAL,
        default=None,
        type=str)
    abfe_check.add_argument(
        '-or', '--out_csv_raw',
        help="The path to output the raw csv file, by default None",
        dest='out_csv_raw',
        nargs=argparse.OPTIONAL,
        default=None,
        type=str)

    args = parser.parse_args()

    if args.command is None:
        print(f"You are using BindFlow: {__version__}")
        print("Chose from the commands options: dag, check_abfe ...")
    elif args.command == "dag":
        dag_maker(input_path=args.input_path, out_name=args.out_name)
    elif args.command == "check_abfe":
        abfe_check_results(
            out_root_folder_path=args.out_root_folder_path,
            out_csv_summary=args.out_csv_summary,
            out_csv_raw=args.out_csv_raw)


if __name__ == "__main__":
    pass
