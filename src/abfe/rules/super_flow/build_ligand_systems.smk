from abfe.preparation import system_builder as sb
import os
import glob
import shutil
from abfe.utils import tools

out_approach_path = config["out_approach_path"]
input_protein_pdb_path = config["inputs"]["protein_pdb_path"]
input_membrane_pdb_path = config["inputs"]["membrane_pdb_path"]

input_ligand_mol_paths = config["inputs"]["ligand_mol_paths"]
ligand_basenames = [os.path.basename(path) for path in input_ligand_mol_paths]
ligand_names = names = [os.path.splitext(ligand_basename)[0] for ligand_basename in ligand_basenames]
# Create a dictionary to map name to basename
ligand_basename_dict = {ligand_name: ligand_basename for ligand_name, ligand_basename in zip(ligand_names, ligand_basenames)}

input_cofactor_mol_path = config["inputs"]["cofactor_mol_path"]
cofactor_on_protein = config["cofactor_on_protein"]

fix_protein = config["fix_protein"]

hmr_factor = float(config['hmr_factor'])

rule make_ligand_copies:
    input:
        input_ligand_mol_paths=input_ligand_mol_paths
    output:
        ligand_copies = expand(out_approach_path+"/{ligand_name}/input/mol/{ligand_basename}", zip, ligand_name=ligand_names, ligand_basename = ligand_basenames)
    run:
        for ligand_name, input_ligand_mol_path, ligand_copy in zip(ligand_names, input.input_ligand_mol_paths, output.ligand_copies):
            shutil.copy(input_ligand_mol_path, ligand_copy)

rule build_ligand_system:
    input:
        mol_file = lambda wildcards: out_approach_path + "/{ligand_name}/input/mol/" + ligand_basename_dict[wildcards.ligand_name]
    output:
        out_approach_path+"/{ligand_name}/input/complex/complex.gro",
        out_approach_path+"/{ligand_name}/input/complex/complex.top",
        # optional(out_approach_path+"/{ligand_name}/input/complex/index.ndx"),
        out_approach_path+"/{ligand_name}/input/ligand/ligand.gro",
        out_approach_path+"/{ligand_name}/input/ligand/ligand.top",
    threads: config["threads"]
    run:
        out_ligand_path = os.path.join(out_approach_path, wildcards.ligand_name)
        out_ligand_input_path = os.path.join(out_ligand_path, 'input')

        # Initialize the files builder
        builder = sb.MakeInputs(
            protein_pdb=input_protein_pdb_path,
            membrane_pdb=input_membrane_pdb_path,
            cofactor_mol=input_cofactor_mol_path,
            cofactor_on_protein=cofactor_on_protein,
            hmr_factor = hmr_factor,
            builder_dir = os.path.join(out_ligand_path, "builder"),
            fix_protein = fix_protein
        )

        # Create topologies and input files
        builder(ligand_mol=input.mol_file,out_dir=out_ligand_input_path)
        builder.clean()