from bindflow.utils import tools
from bindflow.free_energy import mmxbsa_analysis
from pathlib import Path
import tempfile
import os
import shutil


approach_path = config["out_approach_path"]
samples = list(map(str, range(1,1 + config["samples"])))
        

rule run_gmx_mmpbsa:
    input:
        finished = approach_path + "/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/prod.finished",
        top = approach_path + "/{ligand_name}/input/complex/complex.top",
        mmpbsa_in = approach_path + "/{ligand_name}/input/mmpbsa.in",
        ndx = approach_path + "/{ligand_name}/input/complex/index.ndx",
    output:
        mmxbsa_csv = approach_path + "/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/mmxbsa.csv",
    params:
        in_tpr = approach_path + "/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/prod.tpr",
        in_xtc = approach_path + "/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/prod.xtc",
        run_dir = approach_path + "/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/"
    run:
        # Set default host name (reachable as gmx index)
        host_name = 'Protein'
        # Only if debug is activated, update name and selection based on environment variable
        if 'abfe_debug' in os.environ:
            if os.environ['abfe_debug'] == 'True': # environ save the variables as strings 
                if 'abfe_debug_host_name' in os.environ:
                    host_name = os.environ['abfe_debug_host_name']

        # Fix trajectory.
        centered_xtc = tools.center_xtc(
            tpr=params.in_tpr,
            xtc=params.in_xtc,
            run_dir=params.run_dir,
            host_name=host_name
        )

        # Run gmx_mmpbsa
        
        # The index file generated in bindflow.preparation.system_builder.MakeInputs.__call__
        # will always have as first group receptor and as second group ligand
        # therefore, we can pass to the flag -cg <Receptor group> <Ligand group>" = -cg 0 1
        
        gmx_mmpbsa_command = f"gmx_MMPBSA -O -i {input.mmpbsa_in} -cs {params.in_tpr} -ci {input.ndx} -cg 0 1 -ct {params.in_xtc} -cp {input.top} -o res.dat -nogui"
    
        cwd = os.getcwd()
        with tempfile.TemporaryDirectory(prefix='build_', dir=params.run_dir) as tmp_dir:
            try:
                os.chdir(tmp_dir)
                tools.run(gmx_mmpbsa_command)
                mmxbsa_data = mmxbsa_analysis.GmxMmxbsaDataRetriever("COMPACT_MMXSA_RESULTS.mmxsa")
                mmxbsa_data.get_dg().to_csv(output.mmxbsa_csv)
            finally:
                os.chdir(cwd)
                # Clean centered trajectory
                Path(centered_xtc).resolve().unlink()
