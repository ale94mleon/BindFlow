"""Check that BindFlow runs in an environment that has its dependencies.

BindFlow does not declare its dependencies in ``pyproject.toml``: they come from the conda
environment files in ``install/envs``, which the installer uses, and several of them (GROMACS,
AmberTools, the pinned MDRestraintsGenerator, gmx_MMPBSA, the TOFF environment) are not on PyPI.
``pip install bindflow`` therefore installs BindFlow alone. These checks turn a missing dependency
into a message that points to the installer, before anything is submitted.
"""
import importlib.util
import shutil

INSTALL_HINT = (
    "BindFlow's dependencies are installed by its installer, not by pip (pip install bindflow installs "
    "BindFlow alone). See https://bindflow.readthedocs.io/en/latest/source/installation.html"
)


def missing_module_message(module: str) -> str:
    """Return the message for a Python dependency that could not be imported."""
    return f"The Python package '{module}' is missing from this environment. {INSTALL_HINT}"


def check_workflow_dependencies(calculation_type: str, uses_slurm: bool) -> None:
    """Check the dependencies that are only needed once the workflow runs.

    Parameters
    ----------
    calculation_type : str
        ``"fep"`` or ``"mmpbsa"``.
    uses_slurm : bool
        Whether the jobs are submitted to Slurm, which needs the cluster-generic executor plugin.

    Raises
    ------
    RuntimeError
        Listing everything that is missing.
    """
    missing = []
    if shutil.which("snakemake") is None:
        missing.append("the snakemake command")
    modules = {"fep": ["MDRestraintsGenerator"], "mmpbsa": ["GMXMMPBSA"]}.get(calculation_type, [])
    if uses_slurm:
        modules.append("snakemake_executor_plugin_cluster_generic")
    missing += [f"the Python package '{module}'" for module in modules if importlib.util.find_spec(module) is None]
    if missing:
        raise RuntimeError(f"This environment cannot run a {calculation_type} workflow; it is missing "
                           f"{', '.join(missing)}. {INSTALL_HINT}")
