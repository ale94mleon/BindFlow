"""Run small-molecule parameterization (TOFF) and other OpenMM-based tools
from a separate conda environment (by default named ``BindFlow-TOFF``).

The main BindFlow environment does not need OpenMM, OpenFF, RDKit or Espaloma:
those tools are executed out-of-process from the external environment, with a
minimal emulation of ``conda activate``. This works the same whether the
environments were created with micromamba, mamba or conda.

The external environment is located, in this order, from:

#. The ``toff_prefix`` argument (e.g. :func:`bindflow.runners.calculate`).
#. The ``BINDFLOW_TOFF_PREFIX`` environment variable (set by the installer).
#. An environment called ``BindFlow-TOFF`` next to the running one.
"""
import logging
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Sequence, Union

import yaml

logger = logging.getLogger(__name__)

PathLike = Union[os.PathLike, str]

TOFF_ENV_NAME = "BindFlow-TOFF"
TOFF_PREFIX_ENV_VAR = "BINDFLOW_TOFF_PREFIX"
# First TOFF release whose ``parameterize`` CLI accepts ``force_field_type``.
MIN_TOFF_VERSION = "0.3.0"


def _is_toff_env(prefix: Path) -> bool:
    return (prefix / "bin" / "parameterize").is_file()


def find_toff_prefix(toff_prefix: Optional[PathLike] = None) -> Path:
    """Locate the external TOFF environment.

    Parameters
    ----------
    toff_prefix : PathLike, optional
        Explicit prefix of the environment, by default None (auto-discovery)

    Returns
    -------
    Path
        The absolute prefix of the TOFF environment

    Raises
    ------
    RuntimeError
        If the environment cannot be found
    """
    # Explicit locations: if given, they must be valid (no silent fallback)
    explicit = [("toff_prefix argument", toff_prefix),
                (f"${TOFF_PREFIX_ENV_VAR}", os.environ.get(TOFF_PREFIX_ENV_VAR))]
    for source, value in explicit:
        if value:
            prefix = Path(value).expanduser().resolve()
            if _is_toff_env(prefix):
                logger.info(f"Using TOFF environment {prefix} (from {source})")
                return prefix
            raise RuntimeError(f"🚫 {source} points to {prefix}, but {prefix}/bin/parameterize does not exist.")

    # Auto-discovery: an environment named BindFlow-TOFF next to the running one
    main_prefix = Path(sys.prefix).resolve()
    candidates = [main_prefix.parent / TOFF_ENV_NAME,   # <root>/envs/BindFlow -> <root>/envs/BindFlow-TOFF
                  main_prefix / "envs" / TOFF_ENV_NAME]  # running from the base environment
    for prefix in candidates:
        if _is_toff_env(prefix):
            logger.info(f"Using TOFF environment {prefix} (auto-discovered)")
            return prefix

    tried = "\n".join(f"  - {p}" for p in candidates)
    raise RuntimeError(
        f"🚫 Could not find the {TOFF_ENV_NAME} environment. Looked for bin/parameterize in:\n{tried}\n"
        "👉 Run the BindFlow installer, pass toff_prefix=... to bindflow.runners.calculate(), "
        f"or export {TOFF_PREFIX_ENV_VAR}=/path/to/{TOFF_ENV_NAME}.")


def _activated_env(prefix: Path) -> dict:
    """Minimal emulation of ``conda activate`` for executables living in ``prefix``."""
    env = os.environ.copy()
    env["PATH"] = os.pathsep.join([str(prefix / "bin"), env.get("PATH", "")])
    env["CONDA_PREFIX"] = str(prefix)
    env["AMBERHOME"] = str(prefix)   # antechamber/sqm: GAFF and AM1-BCC charges
    env["PYTHONNOUSERSITE"] = "1"    # ignore ~/.local/lib/pythonX.Y/site-packages
    for var in ("PYTHONPATH", "PYTHONHOME"):
        env.pop(var, None)           # would leak the main environment's Python
    return env


def run_in_env(prefix: PathLike, args: Sequence, log_file: Optional[PathLike] = None,
               cwd: Optional[PathLike] = None) -> subprocess.CompletedProcess:
    """Run ``args[0]`` from ``<prefix>/bin`` inside the external environment.

    Parameters
    ----------
    prefix : PathLike
        Prefix of the external environment
    args : Sequence
        The executable name (looked up in ``<prefix>/bin``) followed by its arguments
    log_file : PathLike, optional
        Where to save the command, stdout and stderr, by default None
    cwd : PathLike, optional
        Working directory, by default None

    Returns
    -------
    subprocess.CompletedProcess
        The finished process

    Raises
    ------
    RuntimeError
        If the command fails
    """
    prefix = Path(prefix)
    cmd: List[str] = [str(prefix / "bin" / str(args[0])), *map(str, args[1:])]
    logger.debug(f"Running in {prefix}: {' '.join(cmd)}")
    proc = subprocess.run(cmd, env=_activated_env(prefix), cwd=cwd, capture_output=True, text=True)
    if log_file:
        Path(log_file).write_text(f"$ {' '.join(cmd)}\n\n--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}")
    if proc.returncode != 0:
        tail = "\n".join((proc.stderr or proc.stdout).strip().splitlines()[-20:])
        where = f"\nFull log: {log_file}" if log_file else ""
        raise RuntimeError(f"🚫 Command failed with exit code {proc.returncode}: {' '.join(cmd)}\n{tail}{where}")
    return proc


def _version_tuple(version: str) -> tuple:
    return tuple(int(n) for n in re.findall(r"\d+", version.split("+")[0])[:3])


def check_toff_installation(prefix: PathLike) -> str:
    """Check that TOFF runs in ``prefix`` and is recent enough. Returns its version."""
    proc = run_in_env(prefix, ["parameterize", "--version"])
    version = proc.stdout.strip().split(":")[-1].strip()
    if _version_tuple(version) < _version_tuple(MIN_TOFF_VERSION):
        raise RuntimeError(f"🚫 TOFF {version} found in {prefix}, but BindFlow needs TOFF >= {MIN_TOFF_VERSION}.\n"
                           f"👉 Update it: <conda-frontend> run -p {prefix} python -m pip install -U --no-deps toff")
    logger.info(f"✅ TOFF {version} found in {prefix}")
    return version


def parameterize(prefix: PathLike, input_mol: PathLike, out_dir: PathLike, mol_resi_name: str = "MOL",
                 force_field_type: str = "openff", force_field_code: Optional[str] = None,
                 hmr_factor: Optional[float] = None, safe_naming_prefix: Optional[str] = None,
                 ext_types: Sequence[str] = ("top", "gro")) -> None:
    """Parameterize a small molecule with TOFF from the external environment.

    Writes ``<out_dir>/<mol_resi_name>.<ext>`` for every ``ext`` in ``ext_types``, plus
    ``<mol_resi_name>_toff.yml`` (the TOFF input) and ``<mol_resi_name>_toff.log``.
    The arguments have the same meaning as in :class:`toff.Parameterize`.
    """
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    config = {
        "input_mol": str(Path(input_mol).resolve()),
        "mol_resi_name": mol_resi_name,
        "force_field_type": force_field_type,
        "force_field_code": force_field_code,
        "ext_types": list(ext_types),
        "hmr_factor": float(hmr_factor) if hmr_factor else None,
        "overwrite": True,
        "safe_naming_prefix": safe_naming_prefix,
        "out_dir": str(out_dir),
    }
    config_file = out_dir / f"{mol_resi_name}_toff.yml"
    with open(config_file, "w") as f:
        yaml.safe_dump(config, f)
    run_in_env(prefix, ["parameterize", config_file], log_file=out_dir / f"{mol_resi_name}_toff.log", cwd=out_dir)

    missing = [ext for ext in ext_types if not (out_dir / f"{mol_resi_name}.{ext}").is_file()]
    if missing:
        raise RuntimeError(f"🚫 TOFF finished but did not write {missing} for {input_mol}. "
                           f"See {out_dir / f'{mol_resi_name}_toff.log'}")
