from abfe.preparation import system_builder as sb
import os
import glob
import shutil
from abfe.utils import tools

out_approach_path = config["out_approach_path"]
input_protein_pdb_path = 

ligand_paths = [mol['conf'] for mol in config["inputs"]["ligands"]]
ligand_basenames = [os.path.basename(path) for path in ligand_paths]
ligand_names = names = [os.path.splitext(ligand_basename)[0] for ligand_basename in ligand_basenames]
# Create a dictionary to map name to basename
ligand_dict = {ligand_name: {'basename': ligand_basename, 'definition': ligand_definition} for ligand_name, ligand_basename, ligand_definition in zip(ligand_names, ligand_basenames, config["inputs"]["ligands"])}


hmr_factor = float(config['hmr_factor'])

rule make_ligand_copies:
    input:
        ligand_paths=ligand_paths
    output:
        ligand_copies = expand(out_approach_path+"/{ligand_name}/input/mol/{ligand_basename}", zip, ligand_name=ligand_names, ligand_basename = ligand_basenames)
    run:
        for ligand_path, ligand_copy in zip(input.ligand_paths, output.ligand_copies):
            # TODO: check if the topology was provided and also copy the file
            # I have to use the dict object, I just need this first rule to parallelize the rules
            shutil.copy(ligand_path, ligand_copy)

rule build_ligand_system:
    input:
        # This is just used to paralelize
        mol_file = lambda wildcards: out_approach_path + "/{ligand_name}/input/mol/" + ligand_dict[wildcards.ligand_name]['basename']
    output:
        out_approach_path+"/{ligand_name}/input/complex/complex.gro",
        out_approach_path+"/{ligand_name}/input/complex/complex.top",
        # optional(out_approach_path+"/{ligand_name}/input/complex/index.ndx"),
        out_approach_path+"/{ligand_name}/input/ligand/ligand.gro",
        out_approach_path+"/{ligand_name}/input/ligand/ligand.top",
    params:
        ligand_definition = ligand_dict[wildcards.ligand_name]['definition']
    threads: config["threads"]
    run:
        out_ligand_path = os.path.join(out_approach_path, wildcards.ligand_name)
        out_ligand_input_path = os.path.join(out_ligand_path, 'input')

        # Initialize the files builder
        builder = sb.MakeInputs(
            protein = config["inputs"]["protein"],
            membrane = config["inputs"]["membrane"],
            cofactor = config["inputs"]["cofactor"],
            cofactor_on_protein = config["cofactor_on_protein"],
            water_model = config["water_model"],
            hmr_factor = hmr_factor,
            builder_dir = os.path.join(out_ligand_path, "builder"),
        )

        # Create topologies and input files
        # Here We will use the ligand definition
        builder(ligand_mol=params.ligand_definition,out_dir=out_ligand_input_path)
        builder.clean()