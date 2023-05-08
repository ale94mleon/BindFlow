from abfe.scripts.preparation import system_builder as sb
import os
#Final Check Job
out_approach_path = config["out_approach_path"]
input_protein_pdb_path = config["inputs"]["protein_pdb_path"]
input_membrane_pdb_path = config["inputs"]["membrane_pdb_path"]
input_ligand_mol_paths = config["inputs"]["ligand_mol_paths"]
ligand_names = config['ligand_names']
input_cofactor_mol_path = config["inputs"]["cofactor_mol_path"]
hmr_factor = float(config['hmr_factor'])


rule build_ligand_systems:
    input:
        input_ligand_mol_paths=input_ligand_mol_paths
    output:
        expand(out_approach_path+"/{ligand_name}/input/{sytem_type}/{sytem_type}.{ext}", ligand_name=ligand_names, ext=['gro', 'top'], sytem_type=['ligand', 'complex'])
    run:
        # Initialize the files builder
        builder = sb.MakeInputs(
            protein_pdb=input_protein_pdb_path,
            membrane_pdb=input_membrane_pdb_path,
            cofactor_mol=input_cofactor_mol_path,
            hmr_factor = hmr_factor,
        )

        for input_ligand_mol_path, ligand_name in zip(input['input_ligand_mol_paths'], ligand_names):
            out_ligand_path = os.path.join(out_approach_path, ligand_name)
            out_ligand_input_path = os.path.join(out_ligand_path, 'input')

            # Create topologies and input files
            builder(ligand_mol=input_ligand_mol_path,out_dir=out_ligand_input_path)
        builder.clean()