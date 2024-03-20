import os

root_path = os.path.dirname(__file__)

from bindflow.run_abfe import calculate_abfe
from bindflow.run_mmpbsa import calculate_mmpbsa
from bindflow._version import __version__
