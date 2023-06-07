from abfe.preparation import system_builder as sb
import os
import glob
import shutil
from abfe.utils.tools import makedirs

out_approach_path = config["out_approach_path"]
input_protein_pdb_path = config["inputs"]["protein_pdb_path"]
input_membrane_pdb_path = config["inputs"]["membrane_pdb_path"]
input_ligand_mol_paths = config["inputs"]["ligand_mol_paths"]
ligand_names = config['ligand_names']
input_cofactor_mol_path = config["inputs"]["cofactor_mol_path"]
cofactor_on_protein = config["cofactor_on_protein"]
hmr_factor = float(config['hmr_factor'])

# TODO, build a sanakemake-base parallelizable rule

# It is still not working as expected, if you add new ligands
# it starts the jobs of the already got it results.
rule make_ligand_system_dirs:
    input:
        input_ligand_mol_paths=input_ligand_mol_paths
    output:
        directory(expand(out_approach_path+"/{ligand_name}/input/mol", ligand_name=ligand_names))
    run:
        for ligand_name, input_ligand_mol_path in zip(ligand_names, input.input_ligand_mol_paths):
            makedirs(out_approach_path+f"/{ligand_name}/input/mol/")
            # This is the way that I found to make the next rule parallelizable.
            shutil.copy(input_ligand_mol_path, out_approach_path+f"/{ligand_name}/input/mol/")

rule build_ligand_systems:
    input:
        mol_dir = out_approach_path+"/{ligand_name}/input/mol/"
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
            builder_dir = os.path.join(out_ligand_path, "builder")
        )

        # Create topologies and input files
        new_input_ligand_mol_path = glob.glob(input.mol_dir+"/*")[0]
        builder(ligand_mol=new_input_ligand_mol_path,out_dir=out_ligand_input_path)
        builder.clean()