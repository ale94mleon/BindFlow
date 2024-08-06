import glob
from bindflow.runners import calculate

ligand_mols = glob.glob("inputs/ligands/*mol")

calculate(
    calculation_type='fep',
    protein_pdb_path='inputs/protein.pdb',
    ligand_mol_paths=ligand_mols,
    out_root_folder_path="fep",
    membrane_pdb_path='inputs/membrane.pdb',
    cofactor_mol_path='inputs/dummy_cofactor_23.mol',
    hmr_factor=3,
    approach_name="",
    n_cores_per_job=8,
    num_jobs_receptor_workflow=None,
    num_jobs_per_ligand=40,
    num_replicas=3,
    submit=True,
    cluster_config={})
