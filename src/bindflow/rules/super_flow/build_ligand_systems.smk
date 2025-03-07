from bindflow.preparation import system_builder as sb

ligand_paths = [Path(mol['conf']) for mol in config["inputs"]["ligands"]]
ligand_basenames = [p.name for p in ligand_paths]
ligand_names = [p.stem for p in ligand_paths]
# Create a dictionary to map name to basename
ligand_dict = {ligand_name: {'basename': ligand_basename, 'definition': ligand_definition} for ligand_name, ligand_basename, ligand_definition in zip(ligand_names, ligand_basenames, config["inputs"]["ligands"])}

hmr_factor = config['hmr_factor']
if hmr_factor:
    hmr_factor=float(hmr_factor)
else:
    hmr_factor=None

rule make_ligand_copies:
    input:
        ligand_paths=ligand_paths
    output:
        ligand_copies=expand(out_approach_path+"/{ligand_name}/input/mol/{ligand_basename}", zip, ligand_name=ligand_names, ligand_basename=ligand_basenames)
    run:
        for ligand_path, ligand_copy in zip(input.ligand_paths, output.ligand_copies):
            # TODO: check if the topology was provided and also copy the file
            # I have to use the dict object, I just need this first rule to parallelize the rules
            shutil.copy(ligand_path, ligand_copy)

rule build_ligand_system:
    input:
        # This is just used to parallelize
        mol_file=lambda wildcards: out_approach_path+"/{ligand_name}/input/mol/"+ligand_dict[wildcards.ligand_name]['basename']
    output:
        out_approach_path+"/{ligand_name}/input/complex/complex.gro",
        out_approach_path+"/{ligand_name}/input/complex/complex.top",
        out_approach_path+"/{ligand_name}/input/complex/index.ndx",
        out_approach_path+"/{ligand_name}/input/ligand/ligand.gro",
        out_approach_path+"/{ligand_name}/input/ligand/ligand.top",
    threads: config["threads"]
    run:
        import logging
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(name)s - %(message)s')
        
        out_ligand_path = Path(out_approach_path)/wildcards.ligand_name
        out_ligand_input_path = out_ligand_path/'input'
        
        # Initialize the files builder
        with sb.MakeInputs(
            protein=config["inputs"]["protein"],
            host_name=config["host_name"],
            membrane=config["inputs"]["membrane"],
            cofactor=config["inputs"]["cofactor"],
            cofactor_on_protein=config["cofactor_on_protein"],
            water_model=config["water_model"],
            custom_ff_path=config["custom_ff_path"],
            hmr_factor=hmr_factor,
            fix_protein=config["fix_protein"],
            builder_dir=out_ligand_path/"builder",
            load_dependencies=load_dependencies,
        ) as builder:

            # Create topologies and input files
            # Here We will use the ligand definition
            builder(ligand_definition=ligand_dict[wildcards.ligand_name]['definition'], out_dir=out_ligand_input_path)