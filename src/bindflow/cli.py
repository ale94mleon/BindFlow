#!/usr/bin/env python3

import argparse
import logging
import os

from bindflow import __version__

loggers = [logging.getLogger(name) for name in logging.root.manager.loggerDict]
for logger in loggers:
    logger.setLevel(logging.NOTSET)


def dag_maker(input_path, out_name):
    """Useful for DEBUG"""
    from bindflow.utils import tools
    cwd = os.getcwd()
    os.chdir(input_path)
    tools.run(f"snakemake --dag | dot -Tpng -o {out_name}.png", interactive=True)
    os.chdir(cwd)


def fep_check_results(out_root_folder_path, out_csv_summary, out_csv_raw):
    from bindflow.free_energy import gather_results

    df_summary = gather_results.get_all_fep_dgs(root_folder_path=out_root_folder_path)
    if len(df_summary):
        df_summary = df_summary.sort_values(by='MBAR').reset_index()
        if out_csv_summary:
            df_summary.to_csv(out_csv_summary)
        print(df_summary)
        if out_csv_raw:
            df_raw = gather_results.get_raw_fep_data(root_folder_path=out_root_folder_path)
            if len(df_raw):
                df_raw.to_csv(out_csv_raw)
    else:
        print("🫣")


def mmxbsa_check_results(out_root_folder_path, out_csv_summary, out_csv_raw):
    from bindflow.free_energy import gather_results

    full_df = gather_results.get_raw_mmxbsa_dgs(
        root_folder_path=out_root_folder_path,
        out_csv=out_csv_raw
    )
    df_summary = gather_results.get_all_mmxbsa_dgs(
        full_df=full_df,
        columns_to_process=None,
        out_csv=out_csv_summary
    )
    if len(df_summary):
        print(df_summary)
    else:
        print("🫣")


def clean(out_root_folder_path):
    import subprocess
    import signal
    import shutil
    from pathlib import Path

    print("🧹 Initiating lab cleanup sequence...")

    # ---------- Step 1: Kill Snakemake processes ----------
    try:
        result = subprocess.run(
            ["ps", "aux"], capture_output=True, text=True, check=True
        )
        lines = result.stdout.splitlines()

        # Filter out lines with "snakemake" but not "grep"
        matching = [line for line in lines if "snakemake" in line and "grep" not in line]

        # Extract PIDs (second column)
        pids = [int(line.split()[1]) for line in matching]

        if pids:
            print("🧹🐍 Snakes on the compute node! Initiating containment...")
            for pid in pids:
                os.kill(pid, signal.SIGKILL)
            print("✅ Free of snakes!")
        else:
            print("✅ No Snakemake found slithering around.")
    except subprocess.CalledProcessError:
        print("Failed to run 'ps aux'. Is this a Unix-like system?")
    except IndexError:
        print("Unexpected output format from 'ps'.")
    except ValueError:
        print("Failed to parse PID.")
    except ProcessLookupError:
        print("Process no longer exists.")
    except PermissionError:
        print("Permission denied when trying to kill a process.")
    except Exception as e:
        print("❌ Snake scan failed. Something went wrong:", e)


    # ---------- Step 2: Cancel SLURM jobs ----------
    # Check if 'squeue' is available on the system
    if shutil.which("squeue"):
        print("🧹👨‍🔬 Scanning the SLURM queue...")
        # Build and run the squeue command
        command = 'squeue --noheader -u $USER --format="%i"'
        process = subprocess.Popen(
            command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            executable='/bin/bash'
        )
        output, _ = process.communicate()
        # Decode and filter non-empty job IDs
        job_ids = [jid.strip() for jid in output.decode().split('\n') if jid.strip()]

        if job_ids:
            print("🧹🧪⚠️ Found SLURM jobs. Everything will be cancelled.")
            cancel_command = "scancel " + " ".join(job_ids)
            subprocess.run(cancel_command, shell=True, executable='/bin/bash')
            print("✅ Lab bench cleared! SLURM jobs canceled.")
        else:
            print("✅ No active experiments — lab bench is clear.")
    else:
        print("❌ 'squeue' not found — SLURM job cleanup skipped.")

    # Deleting temporal directories
    # ---------- Step 3: Delete workflow folders ----------

    slurm_logs = Path(out_root_folder_path) / "slurm_logs"
    snakemake =  Path(out_root_folder_path) / ".snakemake"

    if snakemake.exists() and snakemake.is_dir():
        try:
            shutil.rmtree(snakemake)
            print(f"🧽 Removed '{snakemake}' — workflow residue eliminated.")
        except Exception as e:
            print(f"❌ Failed to remove '{snakemake}':", e)
    else:
        print(f"✅ No '{snakemake.name}' directory found — already clean.")

    if slurm_logs.exists() and slurm_logs.is_dir():
        for item in slurm_logs.iterdir():
            item.unlink()
        print(f"🧽 Removed '{slurm_logs}' content  — workflow residue eliminated.")
    else:
        print(f"✅ No '{slurm_logs.name}' directory found — already clean.")

    print("🧼✨ Lab cleanup complete — your workspace is spotless!")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-v', '--version',
        action='version',
        version=f"✨ BindFlow: {__version__}")

    subparsers = parser.add_subparsers(required=True, dest="command")

    dag = subparsers.add_parser('dag',
                                help="🏗️ Build the DAG of the workflow",
                                description="🏗️ Build the DAG of the workflow")
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
    dag.set_defaults(func=lambda args: dag_maker(input_path=args.input_path, out_name=args.out_name))

    fep_check = subparsers.add_parser(
        'check_fep',
        help="🔎 Check for completion of an FEP workflow",
        description="🔎 Check for completion of an FEP workflow")
    fep_check.add_argument(
        dest='out_root_folder_path',
        help='fep directory (`out_root_folder_path` kwarg of :func:`bindflow.runners.calculate`)',
        type=str)
    fep_check.add_argument(
        '-os', '--out_csv_summary',
        help="The path to output the summary csv file, by default None",
        dest='out_csv_summary',
        nargs=argparse.OPTIONAL,
        default=None,
        type=str)
    fep_check.add_argument(
        '-or', '--out_csv_raw',
        help="The path to output the raw csv file, by default None",
        dest='out_csv_raw',
        nargs=argparse.OPTIONAL,
        default=None,
        type=str)
    fep_check.set_defaults(
        func=lambda args: fep_check_results(
            out_root_folder_path=args.out_root_folder_path,
            out_csv_summary=args.out_csv_summary,
            out_csv_raw=args.out_csv_raw))

    mmxbsa_check = subparsers.add_parser(
        'check_mmxbsa',
        help="🔎 Check for completion of an MM(PB/GB)SA workflow",
        description="🔎 Check for completion of an MM(PB/GB)SA workflow")
    mmxbsa_check.add_argument(
        dest='out_root_folder_path',
        help='MM(P/B)BSA directory (`out_root_folder_path` kwarg of :func:`bindflow.runners.calculate`)',
        type=str)
    mmxbsa_check.add_argument(
        '-os', '--out_csv_summary',
        help="The path to output the summary csv file, by default None",
        dest='out_csv_summary',
        nargs=argparse.OPTIONAL,
        default=None,
        type=str)
    mmxbsa_check.add_argument(
        '-or', '--out_csv_raw',
        help="The path to output the raw csv file, by default None",
        dest='out_csv_raw',
        nargs=argparse.OPTIONAL,
        default=None,
        type=str)
    mmxbsa_check.set_defaults(
        func=lambda args: mmxbsa_check_results(
            out_root_folder_path=args.out_root_folder_path,
            out_csv_summary=args.out_csv_summary,
            out_csv_raw=args.out_csv_raw))      

    cleaner = subparsers.add_parser(
        'clean',
        help="🧹 Clean the running directory for restart",
        description="🧹 Clean the running directory for restart. It will:\n"
        "   1 - ⚠️ Kill ALL snakemake process running\n"
        "   2 - ⚠️ Cancel ALL running jobs of an Slurm queue\n"
        "   3 - Delete the .snakemake directory and the content of slurm_logs\n",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    cleaner.add_argument(
        dest='out_root_folder_path',
        help='MM(P/B)BSA / FEP directory (`out_root_folder_path` kwarg of :func:`bindflow.runners.calculate`)',
        type=str)
    cleaner.set_defaults(func=lambda args: clean(out_root_folder_path=args.out_root_folder_path))

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    pass
