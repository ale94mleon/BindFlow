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

rule get_mmxbsa_result:
    input:
        gmxmmpbsa_res = expand(approach_path + "/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/COMPACT_MMXSA_RESULTS.mmxsa", sample = samples)
    output:
        dG_mmxbsa_results = out_approach_path + "/{ligand_names}/{replica}/dG_mmxbsa_results.csv"
    run:
        pass
        # Implementation