#!/usr/bin/env python
import logging
import os
import shutil
import tempfile
from pathlib import Path
from typing import Iterable, List, Tuple, Union

import yaml
from parmed import Structure

from bindflow.home import home
from bindflow.utils import tools

logger = logging.getLogger(__name__)


def get_atom_types(top: tools.PathLike) -> dict:
    """
    Return the atomtypes section as a dict with key atom type and values the
    corresponded line. Include statements are not
    take it into account.
    """
    atom_types = {}

    with open(top, 'r') as f:
        lines = f.readlines()

    section_found = False
    for line in lines:
        if line.startswith('[ atomtypes ]'):
            section_found = True
            continue
        if section_found:
            if line.startswith(';'):
                continue
            elif (not line.strip() or line.startswith('[')) and not line.startswith('[ atomtypes ]'):
                section_found = False
                continue
            fields = line.split()
            if len(fields) >= 6:
                atom_type = fields[0]
                atom_types[atom_type] = line
    return atom_types


def get_molecule_names(input_topology: tools.PathLike, section: str = 'molecules') -> list:
    """It gets the molecule names specified inside input_topology

    Parameters
    ----------
    input_topology : tools.PathLike
        The path of the input topology
    section : str
        The section to extract names from: molecules or moleculetype
    Returns
    -------
    list
        A list of the molecules presented in the topology
    """
    if section not in ['molecules', 'moleculetype']:
        raise ValueError(f"section must be 'molecules', 'moleculetype. {section} was provided.")

    with open(input_topology, 'r') as f:
        lines = f.readlines()

    molecules = []
    i = 0
    while i < len(lines):
        if section in lines[i]:
            i += 1
            while ("[" not in lines[i]):
                if not lines[i].startswith(';'):
                    split_line = lines[i].split()
                    if len(split_line) == 2:
                        molecules.append(split_line[0])
                i += 1
                if i >= len(lines):
                    break
        i += 1

    return molecules


def add_posres_section(input_topology: tools.PathLike, molecules: Iterable[str], out_file: tools.PathLike = None):
    """This will add to the original topology file the corresponded POSRES section to the
    provided molecules:
    Examples of added lines:

    #ifdef POSRES
    #include "posres_{molecule}.itp"
    #endif

    Parameters
    ----------
    input_topology : PathLike
        The path of the input topology
    molecules : Iterable[str]
        The list of name of the molecules for which the topology section will be added
    out_file : PathLike, optional
        The path to output the modified topology file, by default None. Which means that it
        will modify inplace input_topology.
    """
    with open(input_topology, "r") as f:
        top_lines = f.readlines()

    # This is just to be as close as possible to the result of pdb2gmx
    add_sol_internally = False
    if 'SOL' not in molecules:
        molecules.append('SOL')
        add_sol_internally = True

    look_out_flag = False
    out_lines = []
    for line in top_lines:
        if not line.startswith("[ molecules ]"):
            for molecule in molecules:
                if molecule in line and " 3\n" in line:
                    look_out_flag = True
                    mol_name = line.split()[0]
                if look_out_flag and ('[ moleculetype ]' in line or '[ system ]' in line):
                    if mol_name == 'SOL' and add_sol_internally:
                        out_lines.append("\n#ifdef POSRES_WATER\n")
                        out_lines.append("; Position restraint for each water oxygen\n")
                        out_lines.append("[ position_restraints ]\n")
                        out_lines.append(";  i funct       fcx        fcy        fcz\n")
                        out_lines.append("   1    1       1000       1000       1000\n")
                        out_lines.append("#endif\n\n")
                    else:
                        out_lines.append("\n#ifdef POSRES\n")
                        out_lines.append(f'#include "posres_{mol_name}.itp"\n')
                        out_lines.append("#endif\n\n")
                    look_out_flag = False
        out_lines.append(line)

    if not out_file:
        out_file = input_topology

    with open(out_file, "w") as w:
        w.write("".join(out_lines))


def make_posres(input_topology: tools.PathLike, molecules: Iterable[str], out_dir: tools.PathLike, f_xyz: tuple = (2500, 2500, 2500)):
    """Make a position restraint file out of input_topology for all the molecules specified
    on molecules. Taking only the heavy atoms into account

    Parameters
    ----------
    input_topology : PathLike
        The path of the input topology
    molecules : Iterable[str]
        The list of name of the molecules for which the posres file will be created
    out_dir : PathLike
        The path where the posres files will be written
    f_xyz : tuple
        The x, y, z components of the restraint force to be used. It could
        be a float number of a string to be then defined on the mdp file, by default (2500, 2500, 2500)
    """
    for molecule in molecules:
        atom_flag = False

        with open(input_topology, "r") as f:
            top_lines = f.readlines()

        posres_filename = f"posres_{molecule}.itp"
        with open(Path(out_dir)/posres_filename, "w") as posres_file:
            posres_file.write("[ position_restraints ]\n")

            for i in range(len(top_lines)):

                if f"{molecule}  " in (top_lines[i]) and " 3\n" in (top_lines[i]):
                    j = i + 1
                    while j < len(top_lines):

                        if '[ atoms ]' in top_lines[j]:
                            j += 1  # skip this line
                            atom_flag = True

                        if top_lines[j].startswith('['):  # A new section was reached
                            break

                        if atom_flag:
                            if not top_lines[j].startswith("\n") and not top_lines[j].startswith(";") and not top_lines[j].startswith("#"):
                                # Check if heavy atom based on the mass. In case of use of HMR, for that reason 3
                                if float(top_lines[j].split()[7]) > 3:
                                    posres_str = f"{top_lines[j].split()[0]} 1 {f_xyz[0]} {f_xyz[1]} {f_xyz[2]}\n"
                                    posres_file.write(posres_str)
                        j += 1
                    break

    # Add posre sections to the topology
    add_posres_section(input_topology=input_topology, molecules=molecules, out_file=None)


def _tip3p_settles_to_constraints(top: tools.PathLike, molecule: str, out_top: Union[tools.PathLike, None] = None) -> None:
    """Temporal solution to TODO (put the GitHub Issue).
    Basically it will change the settles entrance of `molecule`
    by:
    ; https://gromacs.bioexcel.eu/t/how-to-treat-specific-water-molecules-as-ligand/3470/9
    '[ constraints ]'
    ; ai aj funct length
    1 2 1 0.09572
    1 3 1 0.09572
    2 3 1 0.15139

    Warning
    -------
    This is only useful for replacing the settle section of a tip3p-like molecule.
    This function its just a workaround and will probably bve removed on the future

    Parameters
    ----------
    top : tools.PathLike
        The GMX topology file
    molecule : str
        Name of the molecule where to look for the [ settles ] section
    out_top : Union[tools.PathLike, None], optional
        Path for a output topology, by default None which means that top will be modify in place.
    """
    constraints_section = "; https://gromacs.bioexcel.eu/t/how-to-treat-specific-water-molecules-as-ligand/3470/9\n"\
        "[ constraints ]\n"\
        "; ai aj funct length\n"\
        "1 2 1 0.09572\n"\
        "1 3 1 0.09572\n"\
        "2 3 1 0.15139\n\n"
    with open(top, 'r') as f:
        lines = f.readlines()
    idx_begins, idx_ends = None, None
    section_found = False
    i = 0
    while not lines[i].startswith('[ molecules ]') and i < len(lines):
        if molecule in lines[i] and " 3\n" in lines[i]:
            j = i
            while not lines[j].startswith('[ moleculetype ]') and j < len(lines):
                if lines[j].startswith('[ settles ]'):
                    section_found = True
                    idx_begins = j
                    j += 1
                if section_found and lines[j].startswith(('[', '#')):
                    idx_ends = j
                    break
                j += 1
            break
        i += 1
    if not out_top:
        out_top = top
    with open(out_top, 'w') as f:
        f.write("".join(lines[:idx_begins]) + constraints_section + "".join(lines[idx_ends:]))


class Solvate:
    def __init__(self, water_model_code: str, builder_dir: tools.PathLike = '.solvate', load_dependencies: List[str] = None) -> None:
        """Class to solvate GMX systems.
        Force fields were extracted from `GMX topologies <https://gitlab.com/gromacs/gromacs/-/tree/main/share/top?ref_type=heads>`__.

        Remember to cite properly the main references if you use any of the water models in your work.

        Available water models:
            * amber:
                * amber/spc
                * amber/spce
                * amber/tip3p
                * amber/tip4p
                * amber/tip4pew
                * amber/tip5p
            * charmm
                * charmm/spc
                * charmm/spce
                * charmm/tip3p
                * charmm/tips3p
                * charmm/tip4p
                * charmm/tip5p
            * oplsaa
                * oplsaa/spc
                * oplsaa/spce
                * oplsaa/tip3p
                * oplsaa/tip4p
                * oplsaa/tip4pew
                * oplsaa/tip5p
                * oplsaa/tip5pe

        Parameters
        ----------
        water_model_code : str
            Water model code in the form: "{force field family}/{water model}"
        builder_dir : tools.PathLik, optional
            Where the temporal files will be written.
        load_dependencies : List[str], optional
            It is used in case some previous loading steps are needed for GROMACS commands;
            e.g: ['source /groups/CBG/opt/spack-0.18.1/shared.bash', 'module load sandybridge/gromacs/2022.4'], by default None
        Raises
        ------
        ValueError
            Invalid force field family
        ValueError
            Invalid water model for the selected force field
        """
        self.load_dependencies = load_dependencies

        self.builder_dir = Path(builder_dir).resolve()
        self.builder_dir.mkdir(exist_ok=True, parents=True)

        # Make directory to save topologies after solvation.
        # This will be cleaned out and created every time the class is called (at the beginning)
        # to avoid mismatch.ch between files generated on different calls
        self.solvated_dir = self.builder_dir/'solvated_sys'

        with open(Path(home(dataDir='gmx_water_models'))/'water_models.yml', 'r') as f:
            self.water_models_data = yaml.safe_load(f)

        # Check validity of input code
        force_field_family, water_model = water_model_code.split('/')
        if force_field_family not in self.water_models_data:
            raise ValueError(f"Invalid force field family: {force_field_family}. Choose from {self.water_models_data.keys()}")
        elif water_model not in self.water_models_data[force_field_family]:
            raise ValueError(f"Invalid water model ({water_model}) for {force_field_family}."
                             f"Choose from {self.water_models_data[force_field_family].keys()}")

        self.force_field_family = force_field_family
        self.water_model = water_model
        self.water_itp, self.ions_itp, self.ffnonbonded_itp, self.water_gro = self._get_gmx_water_model()
        self.cwd = os.getcwd()
        if load_dependencies:
            if isinstance(load_dependencies, List):
                self.load_dependencies = tools.HARD_CODE_DEPENDENCIES + ["export GMX_MAXBACKUP=-1"] + load_dependencies
            else:
                raise ValueError(f"load_dependencies must be a List. Provided: {load_dependencies}")

    def _get_gmx_water_model(self) -> Tuple[tools.PathLike]:
        """
        Retrieve water model files

        Returns
        -------
        Tuple[PathLike]
            A tuple with the absolute path of (in this order):
                * water itp file
                * ions itp file
                * water (configuration) gro file
                * atom type itp definition
        """
        # Extract water and ions topologies
        ff_dir = Path(home(dataDir='gmx_water_models'))

        water_itp = (ff_dir/str(self.force_field_family)/f"{self.water_model}.itp").resolve()
        ions_itp = (ff_dir/str(self.force_field_family)/"ions.itp").resolve()
        ffnonbonded_itp = (ff_dir/str(self.force_field_family)/"ffnonbonded.itp").resolve()
        water_gro = (ff_dir/"configurations"/self.water_models_data[self.force_field_family][self.water_model]).resolve()

        return water_itp, ions_itp, ffnonbonded_itp, water_gro

    def _include_all_atom_types(self, top: tools.PathLike) -> None:
        """Add all the atom types of the corresponded force field family
        They will be added to the first [ atomtypes ] section on the file
        top.

        Parameters
        ----------
        top : top: tools.PathLike
            Topology file to be modified.
        """
        with open(top, 'r') as f:
            lines = f.readlines()
        idx_begins, idx_ends = None, None
        section_found = False
        for i, line in enumerate(lines):
            if line.startswith('[ atomtypes ]'):
                idx_begins = i + 1
                section_found = True
                continue
            if section_found:
                if line.startswith('['):
                    idx_ends = i
                    break
        # Add missing atom_types for solvation
        # Extra atom_types will be removed when parmed write the structure
        if idx_begins is not None and idx_ends is not None:
            atom_types = get_atom_types(top)
            for atom_type_name, atom_type_info in get_atom_types(self.ffnonbonded_itp).items():
                if atom_type_name not in atom_types:
                    atom_types[atom_type_name] = atom_type_info
            with open(top, 'w') as f:
                f.write("".join(lines[:idx_begins] + list(atom_types.values()) + ["\n\n"] + lines[idx_ends:]))

    def _include_water_ions_params(self, top: tools.PathLike) -> None:
        """It add include statements to the corresponded water and ion itp files

        Parameters
        ----------
        top : tools.PathLike
            The topology file
        """
        include_statements = [
            f"#include \"{self.water_itp}\"\n",
            f"#include \"{self.ions_itp}\"\n",
        ]
        with open(top, 'r') as f:
            lines = f.readlines()
        for idx, line in enumerate(lines):
            if line.startswith('[ system ]'):
                break
        with open(top, 'w') as f:
            f.write("".join(lines[:idx] + include_statements + lines[idx:]))

    def _add_water_and_ions(
            self,
            gro: tools.PathLike,
            top: tools.PathLike,
            bt: str = "triclinic",
            box: list[float] = None,
            angles: list[float] = None,
            d: float = None,
            c: bool = False,
            pname: str = "NA",
            nname: str = "CL",
            ion_conc: float = 150E-3,
            rmin: float = 1.0) -> None:
        """Make box, solvate and add ions to the system

        Parameters
        ----------
        gro : tools.PathLike
            The configuration file.
        top: tools.PathLike
            The GMX's topology file.
        bt : str, optional
            Box type for -box and -d: triclinic, cubic, dodecahedron, octahedron, by default triclinic
        box : str, optional
            Box vector lengths (a,b,c) in nm (remember that PDB are in Angstroms), by default None. Which means that gmx editconf will use (0 0 0)
        angles : Iterable[float], optional
            This is the angles between the components of the vector in DEGREES.
            It is important that the provided vector has the correct units, by default None.
            For membrane systems (90,90,60) is advisable.
        d : float, optional
            Distance between the solute and the box, by default None. Which means that gmx editconf will use 0
        c : bool, optional
            Center molecule in box (implied by -box and -d), by default False
        pname : str, optional
            Name of the positive ion, by default NA
        nname : str, optional
            Name of the negative ion, by default CL
        ion_conc : float, optional
            Ion concentration used during neutralization of the system, by default 150E-3
        rmin : float, optional
            Minimum distance between ions and non-solvent, by default 1.0
        out_dir : PathLike, optional
            Where the files will be written: solvated.gro, solvated.top, by default '.'
        """
        # We can change directory because all the path used are already converted to absolute paths

        os.chdir(self.solvated_dir)

        editconf_kwargs = dict(
            f=gro,
            o=gro,
            bt=bt
        )
        if box:
            editconf_kwargs['box'] = ' '.join([str(i) for i in box])
        if angles:
            editconf_kwargs['angles'] = ' '.join([str(i) for i in angles])
        if d:
            editconf_kwargs['d'] = d
        if c:
            editconf_kwargs['c'] = True

        # First write an mdp file.
        with open("ions.mdp", "w") as file:
            file.write("; Neighbor searching\n"
                       "cutoff-scheme           = Verlet\n"
                       "rlist                   = 1.1\n"
                       "pbc                     = xyz\n"
                       "verlet-buffer-tolerance = -1\n"
                       "\n; Electrostatics\n"
                       "coulombtype             = cut-off\n"
                       "\n; VdW\n"
                       "rvdw                    = 1.0\n")

        # It is failing becasue There is not define the atom type for the water molecules

        # Define GMX functions
        @tools.gmx_command(load_dependencies=self.load_dependencies)
        def editconf(**kwargs): ...

        @tools.gmx_command(load_dependencies=self.load_dependencies)
        def solvate(**kwargs): ...

        @tools.gmx_command(load_dependencies=self.load_dependencies)
        def grompp(**kwargs): ...

        @tools.gmx_command(load_dependencies=self.load_dependencies, stdin_command="echo \"SOL\"")
        def genion(**kwargs): ...

        # Execute the GMX functions
        editconf(**editconf_kwargs)
        solvate(cp=gro, p=top, cs=self.water_gro, o=gro)
        grompp(f="ions.mdp", c=gro, p=top, o="ions.tpr")
        genion(s="ions.tpr", p=top, o=gro, neutral=True, pname=pname, nname=nname, rmin=rmin, conc=ion_conc)

        # Just to clean the topology. In this way only the used atom types are written.
        # And the include statements are removed
        # It builds a monolithic topology
        struc = tools.readParmEDMolecule(top_file=top, gro_file=gro)
        struc.save(str(top), overwrite=True)
        struc.save(str(gro), overwrite=True)

        # Change back to cwd
        os.chdir(self.cwd)

    def clean(self, directory: Union[None, tools.PathLike] = None) -> None:
        """Used to delete self.builder_dir or directory if provided

        Danger
        ------
        Use it wisely (when directory is provided), you may ended up deleting your computer :-)

        Parameters
        ----------
        directory : Union[None, tools.PathLike], optional
            Directory to delete, by default None
        """
        os.chdir(self.cwd)

        if directory:
            dir2delete = directory
        else:
            dir2delete = self.builder_dir
        try:
            shutil.rmtree(dir2delete)
        except FileNotFoundError:
            pass

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, exception_traceback):
        self.clean()

    def __call__(
            self,
            structure: Structure,
            bt: str = "triclinic",
            box: list[float] = None,
            angles: list[float] = None,
            d: float = None,
            c: bool = False,
            pname: str = "NA",
            nname: str = "CL",
            ion_conc: float = 150E-3,
            rmin: float = 1.0,
            exclusion_list: list = ["SOL", "NA", "CL",],  # "MG", "ZN"],
            out_dir: tools.PathLike = '.',
            out_name: str = 'solvated',
            f_xyz: tuple = (2500, 2500, 2500),
            settles_to_constraints_on: Union[tools.PathLike, str] = None) -> None:

        # Clean any possible generated files during previous calls
        out_dir = Path(out_dir)
        self.clean(self.solvated_dir)
        self.solvated_dir.mkdir(exist_ok=True, parents=True)
        out_dir.mkdir(exist_ok=True, parents=True)

        # Set out files
        init_top = self.solvated_dir/'init.top'
        init_gro = self.solvated_dir/'init.gro'

        # Write the top and gro file of the structure
        structure.save(str(init_top), overwrite=True)
        structure.save(str(init_gro), overwrite=True)

        # Add ion and water section to the topology
        self._include_water_ions_params(init_top)
        # Include all the atom types of the force field
        self._include_all_atom_types(init_top)
        # Solvate add ions and clean the topology
        self._add_water_and_ions(
            gro=init_gro,
            top=init_top,
            bt=bt,
            box=box,
            angles=angles,
            d=d,
            c=c,
            pname=pname,
            nname=nname,
            ion_conc=ion_conc,
            rmin=rmin)
        # Add position restraint section to topology
        molecules = list(set(get_molecule_names(init_top)) - set(exclusion_list))
        # Here I have to move the files to the final directory with its specified names
        make_posres(
            input_topology=init_top,
            molecules=molecules,
            out_dir=out_dir,
            f_xyz=f_xyz
        )
        # Fix conversion for constraints to settles on water-like molecules
        if settles_to_constraints_on:
            _tip3p_settles_to_constraints(top=init_top, molecule=settles_to_constraints_on, out_top=None)

        shutil.copy(init_top, out_dir/f'{out_name}.top')
        shutil.copy(init_gro, out_dir/f'{out_name}.gro')


def index_for_membrane_system(
        configuration_file: tools.PathLike,
        ndxout: tools.PathLike = "index.ndx",
        ligand_name: str = "LIG",
        host_name: str = "Protein",
        cofactor_name: str = None,
        cofactor_on_protein: bool = True,
        load_dependencies: List[str] = None):
    """Make the index file for membrane systems with SOLU, MEMB and SOLV. It uses gmx make_ndx and select internally.
    One examples selection that can be created with ligand_name = LIG; cofactor_name = COF and cofactor_on_protein = True is:
        #. "RECEPTOR" group {host_name};
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
    sele_RECEPTOR = f"\"RECEPTOR\" group {host_name}"
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
    # tools.run(f"""
    #            export GMX_MAXBACKUP=-1
    #            echo "q" | gmx make_ndx -f {configuration_file} -o {tmpndx.name}
    #            gmx select -s {configuration_file} -sf {tmpopt.name} -n {tmpndx.name} -on {ndxout}
    #            """)

    # deleting the line _f0_t0.000 in the file
    with open(ndxout, "r") as index:
        data = index.read()
        data = data.replace("_f0_t0.000", "")
    with open(ndxout, "w") as index:
        index.write(data)

    tmpopt.close()
    tmpndx.close()


def index_for_soluble_system(
        configuration_file: tools.PathLike,
        ndxout: tools.PathLike = "index.ndx",
        ligand_name: str = "LIG",
        host_name: str = "Protein",
        load_dependencies: List[str] = None):
    """Make the index file for soluble system. This is only needed in case MMPBSA calculation;
        #. "RECEPTOR" group {host_name};
        #. "LIGAND" resname {ligand_name};

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
    load_dependencies : List[str], optional
        It is used in case some previous loading steps are needed for GROMACS commands;
        e.g: ['source /groups/CBG/opt/spack-0.18.1/shared.bash', 'module load sandybridge/gromacs/2022.4'], by default None
    """
    tmpopt = tempfile.NamedTemporaryFile(suffix='.opt')
    tmpndx = tempfile.NamedTemporaryFile(suffix='.ndx')

    sele_RECEPTOR = f"\"RECEPTOR\" group {host_name}"
    sele_LIGAND = f"\"LIGAND\" resname {ligand_name}"

    logger.info("Groups in the index.ndx file:")
    logger.info(f"\t{sele_RECEPTOR}")
    logger.info(f"\t{sele_LIGAND}")

    sele_RECEPTOR += ";\n"
    sele_LIGAND += ";\n"

    with open(tmpopt.name, "w") as opt:
        opt.write(sele_RECEPTOR + sele_LIGAND)

    @tools.gmx_command(load_dependencies=load_dependencies, stdin_command="echo \"q\"")
    def make_ndx(**kwargs): ...

    @tools.gmx_command(load_dependencies=load_dependencies)
    def select(**kwargs): ...

    make_ndx(f=configuration_file, o=tmpndx.name)
    select(s=configuration_file, sf=tmpopt.name, n=tmpndx.name, on=ndxout)
    # tools.run(f"""
    #            export GMX_MAXBACKUP=-1
    #            echo "q" | gmx make_ndx -f {configuration_file} -o {tmpndx.name}
    #            gmx select -s {configuration_file} -sf {tmpopt.name} -n {tmpndx.name} -on {ndxout}
    #            """)

    # deleting the line _f0_t0.000 in the file
    with open(ndxout, "r") as index:
        data = index.read()
        data = data.replace("_f0_t0.000", "")
    with open(ndxout, "w") as index:
        index.write(data)

    tmpopt.close()
    tmpndx.close()


if __name__ == '__main__':
    pass
