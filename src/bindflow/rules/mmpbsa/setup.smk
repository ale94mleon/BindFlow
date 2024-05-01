from bindflow.mmpbsa_in.input_loader import GMXMMPBSAInputMaker
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
        if "mmpbsa" in config.keys():
            mmpbsa_config = config["mmpbsa"]
        else:
            mmpbsa_config = {}

        # FIXME: Pass the temperature of the simulation from any mdp file used on the equilibration
        # Ion concentration is used for gb calculation
        # also there is a membrane parameter 
        
        input_file = GMXMMPBSAInputMaker(**mmpbsa_config)
        
        for mmpbsa_in in output.mmpbsa_in:
            input_file.write(mmpbsa_in)

rule mmxbsa_setup:
    input:
        finished = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.finished",
        mdp = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.mdp"
    output:
        mdp = expand(approach_path + "/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/prod.mdp", sample=samples, allow_missing=True),
        gro = expand(approach_path + "/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/init.gro", sample=samples, allow_missing=True)
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
            cmd = f"export GMX_MAXBACKUP=-1; echo \"System\" | gmx trjconv -f {params.in_xtc} -s {params.in_tpr} -o {tmp_dir}/.gro -sep"
            if skip > 0:
                cmd += f" -skip {skip}"
            tools.run(cmd)
            frames = list(Path(tmp_dir).glob('*.gro'))

            if len(frames) >= len(output.gro):
                for frame, gro, mdp_file in zip(frames, output.gro, output.mdp):
                    Path(gro).parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy(frame, gro)
                    mdp_template.write(mdp_file)
            else:
                raise RuntimeError("Not enough frames in equil-mdsim/prod/prod.xtc")