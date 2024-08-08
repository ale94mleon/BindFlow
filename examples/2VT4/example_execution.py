import glob
from bindflow.runners import calculate
from bindflow.orchestration.generate_scheduler import FrontEnd

ligand_mols = glob.glob("inputs/ligands/*mol")
global_config = {
    'cluster': {
        'options': {
            'calculation': None
        }
    }
}

calculate(
    calculation_type='fep',
    protein='inputs/protein.pdb',
    ligands=ligand_mols,
    out_root_folder_path="fep",
    membrane='inputs/membrane.pdb',
    cofactor='inputs/dummy_cofactor_23.mol',
    threads=4,
    num_jobs=12,
    scheduler_class=FrontEnd,
    submit=True,
    global_config=global_config)
