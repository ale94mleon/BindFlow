from pathlib import Path
from bindflow.utils import tools


root_path = Path(__file__).resolve().parent

__PathDir__ = {
    'complex': {
        'membrane': {
            'equi': root_path/'templates/complex/membrane/equi',
            'fep':  root_path/'templates/complex/membrane/fep',
            'mmpbsa': root_path/'templates/complex/membrane/mmpbsa',
        },
        'soluble': {
            'equi': root_path/'templates/complex/soluble/equi',
            'fep':  root_path/'templates/complex/soluble/fep',
            'mmpbsa': root_path/'templates/complex/soluble/mmpbsa',
        },
    },
    'ligand': {
        'equi': root_path/'templates/ligand/equi',
        'fep': root_path/'templates/ligand/fep',
    }
}

_TemplatePath = tools.DotDict(**__PathDir__)
