"""
define restraints for ligand in protein during the uncoupling.
"""

import argparse

import MDAnalysis as mda
from MDRestraintsGenerator import search, restraints

def gen_boresch_restraints():

    parser = argparse.ArgumentParser()
    parser.add_argument('--top', default='../npt_prod1/npt_prod1.tpr',
                        help=('path to input structure topology file '
                                '(e.g. TPR, PRM7)'))
    parser.add_argument('--trj', default='npt_prod1_center.xtc',
                        help=('path to input trajectory file '
                                '(e.g. XTC, NC, TRJ)'))
    parser.add_argument('--ligand_selection', default="resname LIG and not name H*",
                        help='ligand selection string')
    parser.add_argument('--host_selection', default="protein and name CA",
                        help='host atom selection string')
    parser.add_argument('--temperature', type=float, default=298.15,
                        help='simulation temperature')
    parser.add_argument('--outpath', default='./',
                        help='output path for writing files')
    args = parser.parse_args()

    u = mda.Universe(args.top, args.trj)

    # exclude H* named atoms
    ligand_atoms = search.find_ligand_atoms(u, l_selection=args.ligand_selection,
                                            p_align=args.host_selection)

    # find protein atoms
    atom_set = []

    for l_atoms in ligand_atoms:
        psearch = search.FindHostAtoms(u, l_atoms[0],
                                    p_selection=args.host_selection)
        psearch.run(verbose=True)
        atom_set.extend([(l_atoms, p) for p in psearch.host_atoms])

    # Create the boresch finder analysis object
    boresch = restraints.FindBoreschRestraint(u, atom_set)
    boresch.run(verbose=True)

    # boresch.restraint.plot(path=args.outpath) #this is not necessary and might lead to qt errors. (can be turned on if needed)
    boresch.restraint.write(path=args.outpath)

    dG = boresch.restraint.standard_state()

    with open(f'{args.outpath}/dG_off.dat', 'w') as writer:
        writer.write(f'{dG}')


if __name__ == "__main__":...