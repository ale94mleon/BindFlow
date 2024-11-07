from bindflow.free_energy import mmxbsa_analysis
import tempfile

samples = list(map(str, range(1,1 + config["samples"])))
threads = config['threads']
retries = config['retries']
        

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
        import logging
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(name)s - %(message)s')
        logger = logging.getLogger(__name__)
        
        # Fix trajectory.
        centered_xtc = tools.center_xtc(
            tpr=params.in_tpr,
            xtc=params.in_xtc,
            run_dir=params.run_dir,
            host_name=config["host_name"],
            load_dependencies=load_dependencies,
        )

        # Run gmx_mmpbsa
        
        # The index file generated in bindflow.preparation.system_builder.MakeInputs.__call__
        # will always have as first group receptor and as second group ligand
        # therefore, we can pass to the flag -cg <Receptor group> <Ligand group>" = -cg 0 1
        
        frames_for_gmx_mmpbsa_analysis = mdp.get_number_of_frames(params.in_mdp)
        if "mmpbsa" in config.keys():
            mmpbsa_config = config["mmpbsa"]
            if "general" in mmpbsa_config:
                if "startframe" in mmpbsa_config["general"]:
                    if mmpbsa_config["general"]["startframe"] != 0:
                        frames_for_gmx_mmpbsa_analysis = frames_for_gmx_mmpbsa_analysis - mmpbsa_config["general"]["startframe"]
        max_parallel = min(threads, frames_for_gmx_mmpbsa_analysis)
        logger.info(f"📊 Estimated number of frames {frames_for_gmx_mmpbsa_analysis}. Running with {max_parallel} threads.")
        
        dependencies = " && ".join(load_dependencies) + " && "
        gmx_mmpbsa_command = f"gmx_MMPBSA -O -i {input.mmpbsa_in} -cs {params.in_tpr} -ci {input.ndx} -cg 0 1 -ct {centered_xtc} -cp {input.top} -o res.dat -nogui"

        cwd = os.getcwd()
        with tempfile.TemporaryDirectory(prefix='build_', dir=params.run_dir) as tmp_dir:
            os.chdir(tmp_dir)
            try:
                tools.run(f"{dependencies} mpirun --use-hwthread-cpus -np {max_parallel} {gmx_mmpbsa_command}")
            except Exception as e1:
                logger.error(e1)
                try:
                    # Try a second time. We saw that sometimes helps.
                    logger.info(f"🔂 gmx_MMPBSA parallel execution failed; trying a second time...")
                    tools.run(f"{dependencies} mpirun --use-hwthread-cpus -np {max_parallel} {gmx_mmpbsa_command}")
                except Exception as e2:
                    logger.error(e2)
                    logger.info(f"⚠️ gmx_MMPBSA parallel execution failed; switching to sequential execution...")
                    tools.run(dependencies + gmx_mmpbsa_command)
            finally:
                if Path("COMPACT_MMXSA_RESULTS.mmxsa").is_file():
                    logger.info(f"✅ gmx_MMPBSA completed successfully!")
                    mmxbsa_data = mmxbsa_analysis.GmxMmxbsaDataRetriever("COMPACT_MMXSA_RESULTS.mmxsa")
                    mmxbsa_data.store_dg(output.mmxbsa_csv, params.run_dir)
                else:
                    logger.info(f"❌ gmx_MMPBSA execution failed!")

                os.chdir(cwd)
                # Clean centered trajectory
                Path(centered_xtc).resolve().unlink()
