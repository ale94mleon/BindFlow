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


def _is_inside(path, root) -> bool:
    """Check whether ``path`` is ``root`` or lies inside it (both resolved, compared by path parts)."""
    from pathlib import Path
    try:
        Path(path).resolve().relative_to(Path(root).resolve())
        return True
    except (ValueError, OSError):
        return False


def _run_tool(tool: str, *args: str, check: bool = False):
    """Run a system tool by its absolute path and return the completed process.

    The tool is resolved with ``shutil.which`` so the call never depends on a partial path, and the
    arguments are passed as a list without a shell, so nothing in them is interpreted.

    Raises
    ------
    FileNotFoundError
        If ``tool`` is not on the PATH.
    """
    import shutil
    import subprocess  # nosec B404: only used through this helper, with an argument list and no shell

    executable = shutil.which(tool)
    if executable is None:
        raise FileNotFoundError(f"'{tool}' was not found on the PATH")
    # Absolute executable, argument list and no shell; the arguments are PIDs, job IDs and fixed flags.
    return subprocess.run([executable, *args], capture_output=True, text=True, check=check)  # nosec B603


def _process_cwd(pid: int):
    """Return the working directory of a process, or None if it cannot be read (gone, other user, ...)."""
    try:
        return os.readlink(f"/proc/{pid}/cwd")  # Linux
    except FileNotFoundError:
        if not os.path.isdir("/proc/self"):  # No /proc (macOS): ask lsof
            out = _run_tool("lsof", "-a", "-p", str(pid), "-d", "cwd", "-Fn").stdout
            for line in out.splitlines():
                if line.startswith("n"):
                    return line[1:]
    except OSError:
        pass
    return None


def _pipeline_processes(out_root_folder_path) -> list:
    """Return the PIDs of the user's Snakemake processes running inside ``out_root_folder_path``.

    All their descendants (the jobs they started, e.g. gmx mdrun) are included.
    """
    result = _run_tool("ps", "-u", str(os.getuid()), "-o", "pid=,ppid=,args=", check=True)
    children, args = {}, {}
    for line in result.stdout.splitlines():
        parts = line.split(None, 2)
        if len(parts) < 3:
            continue
        pid, ppid = int(parts[0]), int(parts[1])
        args[pid] = parts[2]
        children.setdefault(ppid, []).append(pid)

    own = {os.getpid(), os.getppid()}
    roots = [pid for pid, cmd in args.items()
             if "snakemake" in cmd and pid not in own
             and (cwd := _process_cwd(pid)) and _is_inside(cwd, out_root_folder_path)]

    found, stack = set(), list(roots)
    while stack:
        pid = stack.pop()
        if pid not in found and pid not in own:
            found.add(pid)
            stack.extend(children.get(pid, []))
    return sorted(found)


def _pipeline_slurm_jobs(out_root_folder_path) -> list:
    """Return the IDs of the user's Slurm jobs whose working directory is inside ``out_root_folder_path``.

    BindFlow submits the controller job (job.sh) from that directory and Snakemake submits every
    rule job from there too, so this selects exactly the jobs of this pipeline.
    """
    import getpass
    result = _run_tool("squeue", "--noheader", "-u", getpass.getuser(), "--format=%i|%Z", check=True)
    job_ids = []
    for line in result.stdout.splitlines():
        job_id, _, workdir = line.strip().partition("|")
        if job_id and workdir and _is_inside(workdir, out_root_folder_path):
            job_ids.append(job_id)
    return job_ids


def _stop_pipeline_processes(out_root_folder_path, dry_run: bool) -> None:
    """Stop the Snakemake processes of the pipeline: SIGTERM first, then SIGKILL what is left."""
    import signal
    import time

    try:
        pids = _pipeline_processes(out_root_folder_path)
    except Exception as e:
        print("❌ Snake scan failed. Something went wrong:", e)
        return
    if not pids:
        print("✅ No Snakemake of this pipeline found slithering around.")
        return

    print(f"🧹🐍 {len(pids)} process(es) of this pipeline found: {' '.join(map(str, pids))}")
    if dry_run:
        return
    for sig in (signal.SIGTERM, signal.SIGKILL):
        for pid in pids:
            try:
                os.kill(pid, sig)
            except ProcessLookupError:
                pass
            except PermissionError:
                print(f"Permission denied when trying to stop process {pid}.")
        if sig == signal.SIGTERM:
            time.sleep(3)  # Let Snakemake shut down cleanly
    print("✅ Free of snakes!")


def _cancel_pipeline_slurm_jobs(out_root_folder_path, dry_run: bool) -> None:
    """Cancel the Slurm jobs of the pipeline, if Slurm is available."""
    import shutil

    if not shutil.which("squeue"):
        print("❌ 'squeue' not found — SLURM job cleanup skipped.")
        return
    print("🧹👨‍🔬 Scanning the SLURM queue...")
    try:
        job_ids = _pipeline_slurm_jobs(out_root_folder_path)
    except Exception as e:
        print("❌ Failed to read the SLURM queue:", e)
        return
    if not job_ids:
        print("✅ No active experiments of this pipeline — lab bench is clear.")
        return

    print(f"🧹🧪⚠️ Found {len(job_ids)} SLURM job(s) of this pipeline.")
    if dry_run:
        print("   " + " ".join(job_ids))
        return
    # Chunk to stay below the command-line length limit
    for i in range(0, len(job_ids), 500):
        _run_tool("scancel", *job_ids[i:i + 500])
    print("✅ Lab bench cleared! SLURM jobs canceled.")


def _remove_workflow_folders(out_root_folder_path, dry_run: bool) -> None:
    """Delete the .snakemake directory and the content of slurm_logs."""
    import shutil

    slurm_logs = out_root_folder_path / "slurm_logs"
    snakemake = out_root_folder_path / ".snakemake"
    if dry_run:
        print(f"🔍 Would remove '{snakemake}' and the content of '{slurm_logs}'.")
        return

    if snakemake.is_dir():
        try:
            shutil.rmtree(snakemake)
            print(f"🧽 Removed '{snakemake}' — workflow residue eliminated.")
        except Exception as e:
            print(f"❌ Failed to remove '{snakemake}':", e)
    else:
        print(f"✅ No '{snakemake.name}' directory found — already clean.")

    if slurm_logs.is_dir():
        for item in slurm_logs.iterdir():
            if item.is_dir() and not item.is_symlink():
                shutil.rmtree(item)
            else:
                item.unlink()
        print(f"🧽 Removed '{slurm_logs}' content  — workflow residue eliminated.")
    else:
        print(f"✅ No '{slurm_logs.name}' directory found — already clean.")


def clean(out_root_folder_path, dry_run: bool = False):
    """Clean a pipeline folder for a restart, without touching other pipelines.

    Stops the Snakemake processes running inside the folder (and the jobs they started), cancels the
    Slurm jobs whose working directory is inside it, and deletes the .snakemake directory and the
    content of slurm_logs.

    Parameters
    ----------
    out_root_folder_path : PathLike
        Pipeline folder (``out_root_folder_path`` of :func:`bindflow.runners.calculate`).
    dry_run : bool, optional
        Only list the processes, Slurm jobs and folders that would be affected, by default False.
    """
    from pathlib import Path

    out_root_folder_path = Path(out_root_folder_path).resolve()
    if not out_root_folder_path.is_dir():
        print(f"❌ '{out_root_folder_path}' is not a directory.")
        return

    print(f"🧹 Initiating lab cleanup sequence for '{out_root_folder_path}'...")
    if dry_run:
        print("🔍 Dry run: nothing will be killed, cancelled or deleted.")

    _stop_pipeline_processes(out_root_folder_path, dry_run)
    _cancel_pipeline_slurm_jobs(out_root_folder_path, dry_run)
    _remove_workflow_folders(out_root_folder_path, dry_run)
    if not dry_run:
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
        "   1 - Stop the Snakemake processes running inside the directory (and the jobs they started)\n"
        "   2 - Cancel the Slurm jobs whose working directory is inside the directory\n"
        "   3 - Delete the .snakemake directory and the content of slurm_logs\n"
        "Processes and jobs of other pipelines are not touched. Use --dry-run to see what would be done.\n",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    cleaner.add_argument(
        dest='out_root_folder_path',
        help='MM(P/B)BSA / FEP directory (`out_root_folder_path` kwarg of :func:`bindflow.runners.calculate`)',
        type=str)
    cleaner.add_argument(
        '-n', '--dry-run',
        help="Only list the processes, Slurm jobs and folders that would be affected",
        dest='dry_run',
        action='store_true')
    cleaner.set_defaults(func=lambda args: clean(out_root_folder_path=args.out_root_folder_path, dry_run=args.dry_run))

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    pass
