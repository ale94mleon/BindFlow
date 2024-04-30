import bindflow.mmpbsa_in.mmpbsa_in_load_defaults
from bindflow.utils import tools
from bindflow.mdp import mdp
from bindflow.mdp.templates import TemplatePath
from pathlib import Path
import shutil
import tempfile

# Common to all sub-workflows
approach_path = config["out_approach_path"]
ligand_names = config["ligand_names"]
samples = list(map(str, range(1,1 + config["samples"])))


# TODO
# The parameters should be expose to user customization. A Class like the MDP class may do the job
# Then the new parameters are read from the yml file and update the default values, forcefields and other should be protected and not mutable for MMPBSA (force fields are 
# generated at the beginning.)
rule create_mmxbsa_in:
    output:
        mmpbsa_in = expand(approach_path + "/{ligand_name}/input/mmpbsa.in", ligand_name = ligand_names),
    run:
        import logging
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(name)s - %(message)s')
        
        mmpbsa_config = {}
        if "mmpbsa" in config.keys():
            mmpbsa_config = config["mmpbsa"]
        if not "pb" in mmpbsa_config.keys() and not "gb" in mmpbsa_config.keys():
            mmpbsa_config["pb"] = True

        input_file = bindflow.mmpbsa_in.mmpbsa_in_load_defaults.gmxmmpbsa_input_file(mmpbsa_in_opts = mmpbsa_config)
        
        for mmpbsa_in in output.mmpbsa_in:
            input_file.store_template(mmpbsa_in)

rule mmxbsa_setup:
    input:
        finished = "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.finished",
        mdp = "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.mdp"
    output:
        mdp = expand(approach_path + "/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/prod.mdp", sample=samples),
        gro = expand(approach_path + "/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/init.gro", sample=samples)
    params:
        in_tpr = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.tpr",
        in_xtc = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.xtc",
        sim_dir = approach_path + "/{ligand_name}/{replica}/complex/mmpbsa/simulation",
    run:
        # Update MDP with user provided options
        if config["complex_type"] == 'soluble':
            template_dir = TemplatePath.complex.soluble.mmpbsa
        elif config["complex_type"] == 'membrane':
            template_dir = TemplatePath.complex.membrane.mmpbsa
        mdp_template = mdp.MDP().from_file(template_dir + '/prod.mdp')
        # In case of user defined MDP keywords, take those from the config
        try:
            # TODO sanity check on the passed MDP options
            mdp_template.set_parameters(**config['mdp']['complex']['mmpbsa']['prod'])
        except KeyError:
            pass
        
        sim_dir = Path(params.sim_dir).resolve()
        sim_dir.mkdir(exist_ok=True, parents=True)

        # Get frequency of gro writing
        equil_prod_mdp = mdp.MDP().from_file(input.mdp).parameters
        # The 2 is just to be sure of having enough frames
        skip = int((int(equil_prod_mdp['nsteps'].split(';')[0]) / int(equil_prod_mdp['nstxout-compressed'].split(';')[0])) / (config['samples'] + 2))

        with tempfile.TemporaryDirectory(prefix='split_', dir=sim_dir) as tmp_dir:
            # Get initial configuration from equil/prod
            tools.run(f"export GMX_MAXBACKUP=-1; echo \"System\" | gmx trjconv -f {params.in_xtc} -s {params.in_tpr} -o {tmp_dir}/.gro -sep -skip {skip}")
            frames = list(Path(tmp_dir).glob('*.gro'))

            if len(frames) >= len(output.gro):
                for frame, gro, mdp_file in zip(frames, output.gro, output.mdp):
                    Path(gro).parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy(frame, gro)
                    mdp_template.write(mdp_file)
            else:
                raise RuntimeError("Not enough frames in equil-mdsim/prod/prod.xtc")