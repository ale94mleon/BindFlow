#!/usr/bin/env python3

"""
This script is an example execution running for membrane systems.
"""
import os
import glob
from bindflow import calculate_abfe

ligand_mols = glob.glob("inputs/ligands/*mol")

out_folder = "abfe"

calculate_abfe(
    protein_pdb_path='inputs/protein.pdb',
    ligand_mol_paths=ligand_mols,
    out_root_folder_path="abfe",
    membrane_pdb_path = 'inputs/membrane.pdb',
    cofactor_mol_path = 'inputs/dummy_cofactor_23.mol',
    hmr_factor = 3,
    
    approach_name = "",
    n_cores_per_job= 8,
    num_jobs_receptor_workflow= None,
    num_jobs_per_ligand= 40,
    num_replicas = 3,
    submit= True,
    cluster_config = {})