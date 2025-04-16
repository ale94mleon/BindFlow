import subprocess
import logging


def check_gromacs_installation():
    # Check if the logger is already configured to avoid reconfiguring it
    if not logging.getLogger().hasHandlers():
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    try:
        result = subprocess.run(['gmx', '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode != 0:
            logging.warning("😅 Oops! It seems that GROMACS is not installed or not found in the system PATH. 📦 I really hope that you "
                            "source GROMACS in the global_config[extra_directives][dependencies]. If not, this will get awkward. 🤞")
    except FileNotFoundError:
        logging.warning("😅 Oops! It seems that GROMACS is not installed or not found in the system PATH. 📦 I really hope that you "
                        "source GROMACS in the global_config[extra_directives][dependencies]. If not, this will get awkward. 🤞")
