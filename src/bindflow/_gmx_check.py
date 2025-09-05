import subprocess
import logging
import re


def check_gromacs_installation():
    # Check if the logger is already configured to avoid reconfiguring it
    if not logging.getLogger().hasHandlers():
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    try:
        result = subprocess.run(['gmx', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        if result.returncode == 0:
            installed_version = get_gromacs_version()
            if installed_version is None:
                raise RuntimeError("⚠️ GROMACS was found, but the version could not be determined. "
                                   "Please run `gmx --version` manually and verify your installation.")

            if not is_gromacs_version_leq("2023"):
                raise RuntimeError(
                    f"🚫 Unsupported GROMACS version detected: {installed_version}. "
                    "BinFlow only supports GROMACS versions earlier than 2023 for now."
                    "👉 Please install an older release (e.g., 2021.x or 2022.x) "
                    "or check BindFlow documentation for compatibility details."
                )
            else:
                logging.info(f"✅ Compatible GROMACS version detected: {installed_version}")
        else:
            logging.warning(
                "🤔 Oops! It seems that GROMACS is in the system PATH but failed to run properly. "
                "I really hope that you source GROMACS in the global_config[extra_directives][dependencies]. "
                "If not, this will get awkward. 🤞"
            )
    except FileNotFoundError:
        logging.warning(
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
                logging.warning("⚠️ Could not parse GROMACS version from output.")
                return None
        else:
            return None
    except FileNotFoundError:
        return None


def is_gromacs_version_leq(target_version: str) -> bool:
    """
    Check if the installed GROMACS version is <= target_version.
    target_version should be a string like '2023' or '2023.1'.
    """
    installed_version = get_gromacs_version()
    if installed_version is None:
        return False  # GROMACS not installed or version not detected

    def parse_version(v):
        return tuple(map(int, v.split('.')))

    try:
        return parse_version(installed_version) <= parse_version(target_version)
    except ValueError:
        logging.warning(f"⚠️ Could not compare versions (installed: {installed_version}, target: {target_version})")
        return False