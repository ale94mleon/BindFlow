from ast import parse
import subprocess
import logging
import re


logger = logging.getLogger(__name__)


def check_gromacs_installation():
    try:
        result = subprocess.run(['gmx', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        if result.returncode == 0:
            installed_version = get_gromacs_version()
            if installed_version is None:
                raise RuntimeError("⚠️ GROMACS was found, but the version could not be determined. "
                                   "Please run `gmx --version` manually and verify your installation.")
            min_supported = "2022"
            max_supported = "2026"
            if is_version_in_range(min_version=min_supported, max_version=max_supported):
                logger.info(f"✅ Compatible GROMACS version detected: {installed_version}")
            else:
                raise RuntimeError(
                    f"🚫 Unsupported GROMACS version detected: {installed_version}. "
                    f"Supported versions are >= {min_supported} and < {max_supported}. "
                    "👉 Please install a compatible GROMACS release."
                    )
        else:
            logger.warning(
                "🤔 Oops! It seems that GROMACS is in the system PATH but failed to run properly. "
                "I really hope that you source GROMACS in the global_config[extra_directives][dependencies]. "
                "If not, this will get awkward. 🤞"
            )
    except FileNotFoundError:
        logger.warning(
            "😅 Oops! It seems that GROMACS is not installed or not found in the system PATH. "
            "I really hope that you source GROMACS in the global_config[extra_directives][dependencies]. "
            "If not, this will get awkward. 🤞"
        )


def get_gromacs_version():
    """Return GROMACS version string if installed, otherwise None."""
    try:
        result = subprocess.run(['gmx', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode == 0:
            # Extract version number (usually appears like: "GROMACS version:    2022.3")
            match = re.search(r'GROMACS version:\s*([\d\.]+)', result.stdout)
            if match:
                return match.group(1)
            else:
                logger.warning("⚠️ Could not parse GROMACS version from output.")
                return None
        else:
            return None
    except FileNotFoundError:
        return None


def is_version_in_range(min_version: str | None = None, max_version: str | None = None) -> bool:
    """    Check if the installed GROMACS version lies within
    an optional [min_version, max_version) range.
    target_version should be a string like '2022' or '2022.6'.

    Parameters
    ----------
    min_version : str | None, optional
        inclusive lower bound (>=), by default None
    max_version : str | None, optional
        exclusive upper bound (<), by default None

    Returns
    -------
    bool
        True if installed version is in range
    """
    def parse(v):
        return tuple(map(int, v.split(".")))

    installed_version = get_gromacs_version()
    if installed_version is None:
        return False  # GROMACS not installed or version not detected

    try:
        inst = parse(installed_version)
        if min_version is not None and inst < parse(min_version):
            return False
        if max_version is not None and inst >= parse(max_version):
            return False
        return True
    except ValueError:
        logger.warning(f"⚠️ Could not compare versions (installed: {installed_version}, "
                       f"min={min_version}, max={max_version}")
        return False
