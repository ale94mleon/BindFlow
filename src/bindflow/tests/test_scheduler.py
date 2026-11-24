"""Slurm submission: the generated Snakemake command and job.sh, run against a fake sbatch."""
import importlib.util
import os
import shutil
import stat
import subprocess
import sys

import pytest

from bindflow.orchestration.generate_scheduler import FrontEnd, SlurmScheduler, slurm_validation

CLUSTER = {"partition": "deflt", "time": "01:00:00", "exclusive": True}


def test_slurm_validation_translates_short_keys():
    config = slurm_validation({"p": "deflt", "t": "01:00:00"})
    assert config == {"partition": "deflt", "time": "01:00:00"}


def test_slurm_validation_requires_partition():
    with pytest.raises(ValueError, match="partition"):
        slurm_validation({"time": "01:00:00"})


def test_slurm_validation_rejects_unknown_keys():
    with pytest.raises(ValueError, match="not a valid SLURM"):
        slurm_validation({"partition": "deflt", "not-a-slurm-option": 1})


def test_slurm_command_uses_cluster_generic_executor(tmp_path):
    scheduler = SlurmScheduler(dict(CLUSTER), out_dir=tmp_path, prefix_name="test", snake_executor_file="job.sh")
    command = scheduler.build_snakemake(jobs=10)
    # Snakemake >= 8 removed --cluster and --cluster-cancel
    assert "--executor cluster-generic" in command
    assert "--cluster-generic-cancel-cmd scancel" in command
    assert "--cluster-generic-submit-cmd 'sbatch --parsable" in command
    assert " --cluster " not in command and "--cluster-cancel" not in command
    # Flags are passed bare: sbatch rejects --exclusive=
    assert " --exclusive " in command and "--exclusive=" not in command
    assert "--ntasks={threads}" in command and "--job-name=test.{rule}.{jobid}" in command


def test_slurm_job_script_header(tmp_path):
    scheduler = SlurmScheduler(dict(CLUSTER), out_dir=tmp_path, prefix_name="test", snake_executor_file="job.sh")
    scheduler.build_snakemake(jobs=10)
    scheduler.submit(only_build=True, job_prefix="test")
    lines = (tmp_path / "job.sh").read_text().splitlines()
    assert lines[0] == "#!/bin/bash"
    assert "#SBATCH --partition=deflt" in lines
    assert "#SBATCH --exclusive" in lines
    assert "#SBATCH --job-name=test.RuleThemAll" in lines
    assert lines[-1].startswith("snakemake ")


def test_frontend_command(tmp_path):
    scheduler = FrontEnd(out_dir=tmp_path, snake_executor_file="job.sh")
    command = scheduler.build_snakemake(jobs=4)
    assert command.startswith("snakemake --jobs 4 ")
    assert "cluster" not in command


FAKE_SBATCH = """#!/bin/bash
# Record the arguments, run the job script (the last argument) and print a job ID, as sbatch --parsable does
echo "$@" >> "$SBATCH_LOG"
bash "${@: -1}" > /dev/null 2>&1
echo $RANDOM
"""

SNAKEFILE = """
rule all:
    input: expand("out_{i}.txt", i=[1, 2])

rule make:
    output: "out_{i}.txt"
    threads: 4
    shell: "echo {wildcards.i} > {output}"
"""


@pytest.mark.skipif(sys.platform == "win32", reason="needs bash")
@pytest.mark.skipif(shutil.which("snakemake") is None, reason="needs the snakemake command")
@pytest.mark.skipif(importlib.util.find_spec("snakemake_executor_plugin_cluster_generic") is None,
                    reason="needs snakemake-executor-plugin-cluster-generic")
def test_job_script_runs_with_fake_sbatch(tmp_path):
    """The generated job.sh is a valid Snakemake command: it submits every job through sbatch."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    sbatch = bin_dir / "sbatch"
    sbatch.write_text(FAKE_SBATCH)
    sbatch.chmod(sbatch.stat().st_mode | stat.S_IEXEC)
    workflow = tmp_path / "workflow"
    workflow.mkdir()
    (workflow / "Snakefile").write_text(SNAKEFILE)

    scheduler = SlurmScheduler(dict(CLUSTER), out_dir=workflow, prefix_name="test", snake_executor_file="job.sh")
    scheduler.build_snakemake(jobs=10)
    scheduler.submit(only_build=True, job_prefix="test")

    log = tmp_path / "sbatch.log"
    env = dict(os.environ, PATH=f"{bin_dir}{os.pathsep}{os.environ['PATH']}", SBATCH_LOG=str(log))
    process = subprocess.run([shutil.which("bash"), "job.sh"], cwd=workflow, env=env,  # nosec B603: fixed arguments
                             capture_output=True, text=True, timeout=300)
    assert process.returncode == 0, process.stderr[-2000:]

    assert (workflow / "out_1.txt").read_text().strip() == "1"
    assert (workflow / "out_2.txt").read_text().strip() == "2"
    submissions = log.read_text().splitlines()
    assert len(submissions) == 2
    for submission in submissions:
        assert submission.startswith("--parsable --partition=deflt")
        assert "--ntasks=4" in submission  # {threads} filled in by Snakemake
        assert "--exclusive" in submission.split()
