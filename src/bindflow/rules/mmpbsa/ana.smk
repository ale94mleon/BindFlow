import bindflow.utils.tools
import os
import json
import shutil
import GMXMMPBSA.API
import pandas as pd
import math
import tempfile

approach_path = config["out_approach_path"]
samples = list(map(str, range(1,1 + config["samples"])))
        

rule run_mmxbsa:
    input:
        finished = approach_path + "/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/prod.finished",
        top = approach_path + "/{ligand_name}/input/complex/complex.top",
        mmpbsa_in = approach_path + "/{ligand_name}/input/mmpbsa.in",
        ndx = approach_path + "/{ligand_name}/input/complex/index.ndx",
    output:
        gmxmmpbsa_res = approach_path + "/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/COMPACT_MMXSA_RESULTS.mmxsa",
    params:
        in_tpr = approach_path + "/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/prod.tpr",
        in_xtc = approach_path + "/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/prod.xtc",
        sim_dir = approach_path + "/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/"
    run:
        pass
        # TODO run xtc_center, gmx_mmpbsa and remove XTC center
