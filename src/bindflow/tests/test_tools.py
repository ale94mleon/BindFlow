"""Helpers: tar extraction and archiving, gmx command lines, MDP defaults and topology edits."""
import io
import tarfile
from pathlib import Path

import pytest

from bindflow.mdp import mdp
from bindflow.preparation.solvent import _tip3p_settles_to_constraints
from bindflow.utils import tools


# ---------------------------------------------------------------------------
# Tar extraction and archiving
# ---------------------------------------------------------------------------

def _archive_with(name: str, **attributes) -> tarfile.TarFile:
    """An in-memory tar archive with one member called ``name`` (a file, unless attributes say otherwise)."""
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w") as tar:
        info = tarfile.TarInfo(name)
        data = b"content"
        for key, value in attributes.items():
            setattr(info, key, value)
        info.size = 0 if attributes else len(data)
        tar.addfile(info, None if attributes else io.BytesIO(data))
    buffer.seek(0)
    return tarfile.open(fileobj=buffer)


@pytest.mark.parametrize("name, attributes", [
    ("../escape.txt", {}),
    ("link", {"type": tarfile.SYMTYPE, "linkname": "/etc/passwd"}),
], ids=["path-traversal", "symlink-escape"])
def test_safe_extract_refuses_members_outside_target(tmp_path, name, attributes):
    target = tmp_path / "target"
    # tarfile's "data" filter raises a TarError subclass; the fallback for older Pythons raises ValueError
    with pytest.raises((tarfile.TarError, ValueError)):
        tools.safe_extract_all(_archive_with(name, **attributes), target)
    assert not (tmp_path / "escape.txt").exists()
    assert not (target / "link").exists()


def test_safe_extract_keeps_absolute_paths_inside_target(tmp_path):
    target = tmp_path / "target"
    absolute = tmp_path / "abs.txt"
    archive = _archive_with(str(absolute))
    if hasattr(tarfile, "data_filter"):
        # The "data" filter strips the leading "/" and extracts inside the target
        tools.safe_extract_all(archive, target)
        assert (target / str(absolute).lstrip("/")).read_bytes() == b"content"
    else:
        # The fallback for older Pythons refuses the member
        with pytest.raises(ValueError):
            tools.safe_extract_all(archive, target)
    assert not absolute.exists()


def test_safe_extract_regular_member(tmp_path):
    tools.safe_extract_all(_archive_with("dir/file.txt"), tmp_path)
    assert (tmp_path / "dir" / "file.txt").read_bytes() == b"content"


def test_archive_round_trip_keeps_xtc_inside_project(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    run = Path("project/lig1/1/prod")
    run.mkdir(parents=True)
    (run / "prod.xtc").write_text("trajectory")
    (run / "prod.mdp").write_text("nsteps = 10")
    tools.archive("project", name="archive", out_check_file=False)
    with tarfile.open("archive.tar") as tar:
        assert "lig1/1/prod/prod.xtc" in tar.getnames()  # relative to the project, not an absolute path
    tools.unarchive("archive.tar", "restored")
    assert (tmp_path / "restored/lig1/1/prod/prod.xtc").read_text() == "trajectory"
    assert (tmp_path / "restored/lig1/1/prod/prod.mdp").read_text() == "nsteps = 10"


# ---------------------------------------------------------------------------
# gmx command lines
# ---------------------------------------------------------------------------

def test_gmx_command_renders_options(monkeypatch):
    commands = []
    monkeypatch.setattr(tools, "run", lambda command, **kwargs: commands.append(command))

    @tools.gmx_command(load_dependencies=["module load gromacs"])
    def mdrun(**kwargs): ...

    mdrun(deffnm="prod", gpu_id=0, cpi=True, nb=None, rerun=False, pin="")
    # 0 is a value (gpu_id: 0 used to be dropped); None, False and "" mean "not set"; True is a bare flag
    assert commands == ["module load gromacs && gmx mdrun -deffnm prod -gpu_id 0 -cpi"]


# ---------------------------------------------------------------------------
# MDP defaults
# ---------------------------------------------------------------------------

def test_mdp_parameters_do_not_leak_between_instances():
    default_nsteps = mdp.MDP().parameters.get("nsteps")
    mdp.MDP(nsteps=5)
    assert mdp.MDP().parameters.get("nsteps") == default_nsteps
    assert mdp._MDP_PARAM_DEFAULT.get("nsteps") == default_nsteps


# ---------------------------------------------------------------------------
# [ settles ] -> [ constraints ]
# ---------------------------------------------------------------------------

WATER = "[ moleculetype ]\nSOL  3\n\n"
SETTLES = "[ settles ]\n1 1 0.09572 0.15139\n\n"


@pytest.mark.parametrize("after_settles", ["[ exclusions ]\n1 2 3\n\n[ molecules ]\nSOL 10\n", ""],
                         ids=["section-follows", "end-of-file"])
def test_settles_replaced_by_constraints(tmp_path, after_settles):
    top = tmp_path / "topol.top"
    top.write_text(WATER + SETTLES + after_settles)
    _tip3p_settles_to_constraints(top, "SOL")
    text = top.read_text()
    assert "[ settles ]" not in text
    assert text.count("[ constraints ]") == 1
    assert text.count("[ moleculetype ]") == 1  # the topology is not written twice


@pytest.mark.parametrize("content", [WATER + "[ atoms ]\n1 OW\n\n[ molecules ]\nSOL 10\n", WATER],
                         ids=["no-settles", "no-molecules-section"])
def test_settles_missing_raises_and_keeps_file(tmp_path, content):
    top = tmp_path / "topol.top"
    top.write_text(content)
    with pytest.raises(ValueError, match="settles"):
        _tip3p_settles_to_constraints(top, "SOL")
    assert top.read_text() == content
