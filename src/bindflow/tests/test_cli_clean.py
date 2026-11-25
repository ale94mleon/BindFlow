"""bindflow clean: only the processes and Slurm jobs of the given pipeline are touched."""
import os
import shutil
import signal
import stat
import subprocess  # nosec B404: starts fake snakemake processes for clean to find
import sys
import time

import pytest

from bindflow.cli import clean

pytestmark = pytest.mark.skipif(sys.platform == "win32", reason="needs a POSIX shell")


def _executable(path, content):
    path.write_text(content)
    path.chmod(path.stat().st_mode | stat.S_IEXEC)


@pytest.fixture()
def pipelines(tmp_path, monkeypatch):
    """Two pipelines and a folder whose name starts like the first, with fake squeue and scancel on the PATH."""
    for name in ("pipeA", "pipeB", "pipeAB"):
        (tmp_path / name).mkdir()
    (tmp_path / "pipeA" / ".snakemake").mkdir()
    (tmp_path / "pipeA" / "slurm_logs" / "sub").mkdir(parents=True)
    (tmp_path / "pipeA" / "slurm_logs" / "job.out").write_text("log")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    # Job 102 runs in a subfolder of pipeA, as rule jobs do
    _executable(bin_dir / "squeue", f"""#!/bin/sh
echo "101|{tmp_path}/pipeA"
echo "102|{tmp_path}/pipeA/lig1/1"
echo "103|{tmp_path}/pipeB"
echo "104|{tmp_path}/pipeAB"
""")
    _executable(bin_dir / "scancel", f'#!/bin/sh\necho "$@" >> "{tmp_path}/scancel.log"\n')
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    return tmp_path


def test_clean_dry_run_touches_nothing(pipelines, capsys):
    clean(pipelines / "pipeA", dry_run=True)
    out = capsys.readouterr().out
    assert "101 102" in out
    assert not (pipelines / "scancel.log").exists()
    assert (pipelines / "pipeA" / ".snakemake").is_dir()
    assert (pipelines / "pipeA" / "slurm_logs" / "job.out").exists()


def test_clean_cancels_only_the_pipeline_jobs(pipelines):
    clean(pipelines / "pipeA")
    assert (pipelines / "scancel.log").read_text().split() == ["101", "102"]
    assert not (pipelines / "pipeA" / ".snakemake").exists()
    assert list((pipelines / "pipeA" / "slurm_logs").iterdir()) == []  # subfolders included


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="the fake process name relies on Linux ps")
def test_clean_stops_only_the_pipeline_processes(pipelines):
    def fake_snakemake(folder):
        # A process whose command line contains "snakemake", running in folder, with a child process
        # (own process group, so the test can clean up the child too)
        command = [shutil.which("bash"), "-c", "exec -a 'snakemake --jobs 1' bash -c 'sleep 60 & wait'"]
        return subprocess.Popen(command, cwd=folder, start_new_session=True)  # nosec B603: fixed arguments

    inside, outside = fake_snakemake(pipelines / "pipeA"), fake_snakemake(pipelines / "pipeB")
    try:
        time.sleep(0.5)
        clean(pipelines / "pipeA")
        assert inside.wait(timeout=10) is not None
        assert outside.poll() is None, "a process of another pipeline was stopped"
    finally:
        for process in (inside, outside):
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()
