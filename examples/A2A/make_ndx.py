from pathlib import Path
from glob import glob
from bindflow.utils import tools
from typing import List
import tempfile
import logging
logger = logging.getLogger(__name__)


def index_for_membrane_system(
        configuration_file: tools.PathLike,
        ndxout: tools.PathLike = "index.ndx",
        ligand_name: str = "LIG",
        host_name: str = "Protein",
        host_name_for_mmpbsa: str = "Protein",
        cofactor_name: str = None,
        cofactor_on_protein: bool = True,
        load_dependencies: List[str] = None):
    """Make the index file for membrane systems with SOLU, MEMB and SOLV. It uses gmx make_ndx and select internally.
    One examples selection that can be created with ligand_name = LIG; cofactor_name = COF and cofactor_on_protein = True is:
        #. "RECEPTOR" group {host_name_for_mmpbsa};
        #. "LIGAND" resname {ligand_name};
        #. "SOLU" group {host_name} or resname {ligand_name} or resname COF;
        #. "MEMB" ((group System and ! group Water_and_ions) and ! group {host_name}) and ! (resname {ligand_name}) and ! (resname COF);
        #. "SOLV" group Water_and_ions;


    Parameters
    ----------
    configuration_file : PathLike
        PDB or GRO file of the system.
    ndxout : PathLike
        Path to output the index file.
    ligand_name : str
        The residue name for the ligand in the configuration file, by default "LIG".
    host_name : str
        The group name for the host in the configuration file, by default "Protein".
    cofactor_name : str
        The residue name for the cofactor in the configuration file, bt default None
    cofactor_on_protein : bool
        It only works if cofactor_name is provided. If True, the cofactor will be part of the protein and the lignad
        if False will be part of the solvent and ions, bt default True
    load_dependencies : List[str], optional
        It is used in case some previous loading steps are needed for GROMACS commands;
        e.g: ['source /groups/CBG/opt/spack-0.18.1/shared.bash', 'module load sandybridge/gromacs/2022.4'], by default None
    """
    tmpopt = tempfile.NamedTemporaryFile(suffix='.opt')
    tmpndx = tempfile.NamedTemporaryFile(suffix='.ndx')
    # Nice use of gmx select, see the use of the parenthesis
    sele_RECEPTOR = f"\"RECEPTOR\" group {host_name_for_mmpbsa}"
    sele_LIGAND = f"\"LIGAND\" resname {ligand_name}"
    sele_MEMB = f"\"MEMB\" ((group System and ! group Water_and_ions) and ! group {host_name}) and ! (resname {ligand_name})"
    sele_SOLU = f"\"SOLU\" group {host_name} or resname {ligand_name}"
    sele_SOLV = "\"SOLV\" group Water_and_ions"
    if cofactor_name:
        sele_MEMB += f" and ! (resname {cofactor_name})"
        if cofactor_on_protein:
            sele_SOLU += f" or resname {cofactor_name}"
        else:
            sele_SOLV += f" or resname {cofactor_name}"

    logger.info("Groups in the index.ndx file:")
    logger.info(f"\t{sele_RECEPTOR}")
    logger.info(f"\t{sele_LIGAND}")
    logger.info(f"\t{sele_SOLU}")
    logger.info(f"\t{sele_MEMB}")
    logger.info(f"\t{sele_SOLV}")

    sele_RECEPTOR += ";\n"
    sele_LIGAND += ";\n"
    sele_SOLU += ";\n"
    sele_MEMB += ";\n"
    sele_SOLV += ";\n"

    with open(tmpopt.name, "w") as opt:
        opt.write(sele_RECEPTOR + sele_LIGAND + sele_SOLU + sele_MEMB + sele_SOLV)

    @tools.gmx_command(load_dependencies=load_dependencies, stdin_command="echo \"q\"")
    def make_ndx(**kwargs): ...

    @tools.gmx_command(load_dependencies=load_dependencies)
    def select(**kwargs): ...

    make_ndx(f=configuration_file, o=tmpndx.name)
    select(s=configuration_file, sf=tmpopt.name, n=tmpndx.name, on=ndxout)

    # deleting the line _f0_t0.000 in the file
    with open(ndxout, "r") as index:
        data = index.read()
        data = data.replace("_f0_t0.000", "")
    with open(ndxout, "w") as index:
        index.write(data)

    tmpopt.close()
    tmpndx.close()


index_files = glob('mmpbsa/*/*/input/complex/index.ndx')

for index_file in index_files:
    index_file = Path(index_file)
    parent_dir = index_file.parent
    print(index_file)

    index_for_membrane_system(
        configuration_file=parent_dir/'complex.gro',
        ndxout=index_file,
        ligand_name='LIG',
        host_name='Protein',
        host_name_for_mmpbsa='Protein or resname COF',
        load_dependencies=['export GMX_MAXBACKUP=-1'],
        cofactor_name='COF',
        cofactor_on_protein=True,
    )
