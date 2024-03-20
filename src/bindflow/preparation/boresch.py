"""
define restraints for ligand in protein during the uncoupling.
"""

import MDAnalysis as mda
from MDRestraintsGenerator import restraints, search

from bindflow.utils.tools import PathLike


def gen_restraint(topology: PathLike,
                  trajectory: PathLike,
                  ligand_selection: str = 'resname LIG and not name H*',
                  host_selection: str = 'protein and name CA',
                  temperature: float = 298.15,
                  outpath: PathLike = './'):
    """It will generate the Boresch restraints. It use MDAnalysis and MDRestraintsGenerator.
    It defines restraints for ligand in protein during the uncoupling.

    Parameters
    ----------
    topology : PathLike
        Path to the topology (binary) file. E.g: TPR, PRM7
    trajectory : PathLike
        Path to the trajectory file. E.g: XTC, NC, TRJ
    ligand_selection : str, optional
        MDAnalysis selection to define the ligand, by default 'resname LIG and not name H*'
    host_selection : str, optional
        MDAnalysis selection to define the host (receptor), by default 'protein and name CA'
    temperature : float
        simulation temperature [298.15]
    outpath : PathLike, optional
        Where the output files will be written out, by default './'
    """

    u = mda.Universe(topology, trajectory)

    # exclude H* named atoms
    ligand_atoms = search.find_ligand_atoms(u, l_selection=ligand_selection,
                                            p_align=host_selection)

    # find protein atoms
    atom_set = []

    for l_atoms in ligand_atoms:
        psearch = search.FindHostAtoms(u, l_atoms[0], p_selection=host_selection)
        psearch.run(verbose=True)
        atom_set.extend([(l_atoms, p) for p in psearch.host_atoms])

    # Create the boresch finder analysis object
    boresch = restraints.FindBoreschRestraint(u, atom_set)
    boresch.run(verbose=True)

    # boresch.restraint.plot(path=args.outpath) #this is not necessary and might lead to qt errors. (can be turned on if needed)
    boresch.restraint.write(path=outpath)

    dG = boresch.restraint.standard_state(temperature=temperature)

    with open(f'{outpath}/dG_off.dat', 'w') as writer:
        writer.write(f'{dG}')


if __name__ == "__main__":
    pass
