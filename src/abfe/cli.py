#!/usr/bin/env python3

import glob, os
import argparse

from abfe import calculate_abfe, __version__
from abfe.utils import tools
import logging


loggers = [logging.getLogger(name) for name in logging.root.manager.loggerDict]
for logger in loggers:
    logger.setLevel(logging.NOTSET)

# TODO, fix this cli option is NOT UPDATE
def abfe_run():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', "--protein_pdb_path", help='Input protein pdb file path', required=True, type=str)
    parser.add_argument('-l', "--ligand_mol_dir",  help='Input ligand(s) mol file path', required=True, type=str) 
    parser.add_argument('-o', "--output_dir_path",  help='Output approach folder', required=True, type=str)
    parser.add_argument('-c', "--cofactor_mol_path", help='Input cofactor(s) mol file path', required=False, default=None, type=str)
    parser.add_argument('-m', "--membrane_pdb_path",
                        help="Input membrane pdb file path. WARNING: The CRYST1 information of this PDB will be used for solvating the system."\
                            "The protein-membrane system MUST be aligned to the Z-axis (needed for pressure couple scheme)."\
                            "CHARMM-GUI is a good option to get this file.",
                        required=False, default=None, type=str)
    parser.add_argument('-hmr', "--hmr_factor",
                        help='The Hydrogen Mass Repartition factor to use. 4 fs of integration time step will be used no matter what hmf_factor is provided. '\
                            'Values greater than 2 are advised, if not the system may be unstable.',
                        required=False, default=3.0, type=float)
    parser.add_argument('-t', "--threads", help='Maximum number of threads to us on the rules', required=False, default=12, type=int)
    parser.add_argument('-lj', "--ligand_jobs", help='Number of jobs in parallel for receptor workflow. By defaults it will take number of ligands * replicas', required=False, default=None, type=int)
    parser.add_argument('-jlj', "--jobs_per_ligand_job", help='Number of jobs in parallel for each ligand workflow. By defaults it will take number of ligands * replicas', required=False, default=10000, type=int)
    parser.add_argument('-r', "--replicas", help='Number of replicates', required=False, default=3, type=int)
    parser.add_argument('-gc', '--global_config',
                        help="This is the configuration YAML file for your simulation. Here you define the cluster characteristics and other optional parameters of the simulation", required=True, type=str)
    parser.add_argument('-submit',help='Will automatically submit the ABFE calculations', required=False, action='store_true')
    parser.add_argument('-v', '--version', action='version', version=f"abfe: {__version__}")

    args = parser.parse_args()
    print(args)
    mol_paths = [f for pattern in ["*.mol", "*.sdf"] for f in glob.glob(os.path.join(args.ligand_mol_dir, pattern))]
    
    calculate_abfe(protein_pdb_path = args.protein_pdb_path,
                   ligand_mol_paths = mol_paths,
                   out_root_folder_path=args.output_dir_path,
                   cofactor_mol_path = args.cofactor_mol_path,
                   membrane_pdb_path = args.membrane_pdb_path,
                   hmr_factor = args.hmr_factor,
                   threads = args.threads,
                   ligand_jobs = args.ligand_jobs,
                   jobs_per_ligand_job =args.jobs_per_ligand_job,
                   replicas = args.replicas,
                   submit = args.submit,
                   global_config=args.global_config,
                   )

def abfe_dag():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', "--path", help = 'Input protein pdb file path', type=str, default = '.')
    parser.add_argument('-o', "--out_name",  help = 'Input ligand(s) mol file path', default = 'dag', type = str) 
    
    args = parser.parse_args()
    print(args)
    cwd = os.getcwd()
    os.chdir(args.path)
    tools.run(f"snakemake --dag | dot -Tpng -o {args.out_name}.png", interactive=True)
    os.chdir(cwd)

def abfe_check_results():
    from abfe.free_energy import gather_results
    parser = argparse.ArgumentParser()
    parser.add_argument(
        dest = 'root_folder_path',
        help='Input protein pdb file path',
        type=str)
    args = parser.parse_args()

    print(gather_results.get_all_dgs(root_folder_path=args.root_folder_path))
if __name__ == "__main__":...
