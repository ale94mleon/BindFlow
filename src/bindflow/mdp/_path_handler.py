import os

from bindflow.utils import tools

root_path = os.path.abspath(os.path.dirname(__file__))

__PathDir__ = {
    'complex': {
        'membrane': {
            'equi': os.path.join(root_path, 'templates', 'complex', 'membrane', 'equi'),
            'fep':  os.path.join(root_path, 'templates', 'complex', 'membrane', 'fep'),
            'mmpbsa': os.path.join(root_path, 'templates', 'complex', 'membrane', 'mmpbsa'),
        },
        'soluble': {
            'equi': os.path.join(root_path, 'templates', 'complex', 'soluble', 'equi'),
            'fep':  os.path.join(root_path, 'templates', 'complex', 'soluble', 'fep'),
            'mmpbsa': os.path.join(root_path, 'templates', 'complex', 'soluble', 'mmpbsa'),
        },
    },
    'ligand': {
        'equi': os.path.join(root_path, 'templates', 'ligand', 'equi'),
        'fep': os.path.join(root_path, 'templates', 'ligand', 'fep'),
    }
}

_TemplatePath = tools.DotDict(**__PathDir__)
