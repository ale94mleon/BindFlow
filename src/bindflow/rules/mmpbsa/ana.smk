from bindflow.free_energy import mmxbsa_analysis
import tempfile

samples = list(map(str, range(1,1 + config["samples"])))
threads = config['threads']
        

rule run_gmx_mmpbsa:
    input:
        finished=out_approach_path+"/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/prod.finished",
        top=out_approach_path+"/{ligand_name}/input/complex/complex.top",
        mmpbsa_in=out_approach_path+"/{ligand_name}/input/mmpbsa.in",
        ndx=out_approach_path+"/{ligand_name}/input/complex/index.ndx",
    output:
        mmxbsa_csv=out_approach_path+"/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/mmxbsa.csv",
    params:
        in_tpr=out_approach_path+"/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/prod.tpr",
        in_xtc=out_approach_path+"/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/prod.xtc",
        in_mdp=out_approach_path+"/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/prod.mdp",
        run_dir=out_approach_path+"/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/"
    threads: threads
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
        
        max_parallel = min(threads, mdp.get_number_of_frames(params.in_mdp))
        print(f"Estimated number of frames {mdp.get_number_of_frames(params.in_mdp)} are run with {max_parallel} threads.")
        gmx_mmpbsa_command = f"mpirun -np {max_parallel} gmx_MMPBSA -O -i {input.mmpbsa_in} -cs {params.in_tpr} -ci {input.ndx} -cg 0 1 -ct {centered_xtc} -cp {input.top} -o res.dat -nogui"
    
        cwd = os.getcwd()
        with tempfile.TemporaryDirectory(prefix='build_', dir=params.run_dir) as tmp_dir:
            try:
                os.chdir(tmp_dir)
                tools.run(gmx_mmpbsa_command)
                mmxbsa_data = mmxbsa_analysis.GmxMmxbsaDataRetriever("COMPACT_MMXSA_RESULTS.mmxsa")
                mmxbsa_data.store_dg(output.mmxbsa_csv, params.run_dir)
            finally:
                os.chdir(cwd)
                # Clean centered trajectory
                Path(centered_xtc).resolve().unlink()
