#!/user/bin python
import os
import tempfile
from abfe.utils import tools
from abfe.utils.tools import PathLike
from typing import Iterable

def get_molecule_names(input_topology:PathLike, section:str = 'molecules') -> list:
    """It gets the molecule names specified inside input_topology

    Parameters
    ----------
    input_topology : PathLike
        The path of the input topology
    input_topology : str
        The section to extract names from molecules or moleculetype
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
                if i >= len(lines): break
        i += 1

    return molecules

def add_posres_section(input_topology:PathLike, molecules:Iterable[str], out_file:PathLike="topol2.top"):
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
        The path to output the modified topology file, by default "topol2.top"
    """
    with open(input_topology, "r") as f:
        top_lines = f.readlines()

    look_out_flag = False
    out_line = []
    for line in top_lines:
        for molecule in molecules:
            if f"{molecule}" in line and " 3\n" in line:
                look_out_flag = True
                mol_name = line.split()[0]
            if look_out_flag and ('[ moleculetype ]' in line or '[ system ]' in line):
                out_line.append("\n#ifdef POSRES\n")
                out_line.append(f'#include "posres_{mol_name}.itp"\n')
                out_line.append("#endif\n\n")
                # if ("LIG" in mol_name or "MOL" in mol_name):
                #     out_line.append("\n#ifdef POSRES_LIG\n")
                #     out_line.append(f'#include "posres_{mol_name}.itp"\n')
                #     out_line.append("#endif\n\n")
                look_out_flag = False
        out_line.append(line)

    with open(out_file, "w") as w:
        for line in out_line:
            w.write(line)

def make_posres_files(input_topology:PathLike, molecules:Iterable[str], out_dir:PathLike, f_xyz:tuple = (2500, 2500, 2500)):
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
        bonds_flag = False

        with open(input_topology, "r") as f:
            top_lines = f.readlines()

        posres_filename = f"posres_{molecule}.itp"
        with open(os.path.join(out_dir,posres_filename), "w") as posres_file:
            posres_file.write("[ position_restraints ]\n")

            for i in range(len(top_lines)):
                if f"{molecule}  " in (top_lines[i]) and " 3\n" in (top_lines[i]):
                    atom_flag = False
                    bonds_flag = False

                    for j in range(i + 1, len(top_lines)):
                        if '[ atoms ]' in top_lines[j]:
                            atom_flag = True
                        if '[ bonds ]' in top_lines[j]:
                            bonds_flag = True
                            break
                        if atom_flag and not bonds_flag:
                            if not "[" in top_lines[j] and not top_lines[j].startswith("\n") and not top_lines[j].startswith(";"):
                                # Check if heavy atom based on the mass. In case of use of HMR, for that reason 3
                                if float(top_lines[j].split()[7]) > 3:
                                    posres_str = f"{top_lines[j].split()[0]} 1 {f_xyz[0]} {f_xyz[1]} {f_xyz[2]}\n"
                                    posres_file.write(posres_str)

def make_ion_moleculetype_section(ion_name:str) -> str:
    """A simple function to create the moleculetype section
    for an ion

    Parameters
    ----------
    ion_name : str
        The ion name, only valid: CL, K, NA

    Returns
    -------
    str
        The GROMACS topology [ moleculetype ] section

    Raises
    ------
    ValueError
        If invalid ion_name
    """
    internal_data = {
        'CL':{
            'NAME':'CL',
            'TYPE':'Cl',
            'RESIDUE':'CL',
            'ATOM':'CL',
            'CHARGE':-1.000000,
            'MASS':35.450000
        },
        'K':{
            'NAME':'CL',
            'TYPE':'Cl',
            'RESIDUE':'CL',
            'ATOM':'CL',
            'CHARGE':1.000000,
            'MASS':39.100000
        },
        'NA':{
            'NAME':'NA',
            'TYPE':'Na',
            'RESIDUE':'NA',
            'ATOM':'NA',
            'CHARGE':1.000000,
            'MASS':22.990000
        }, 
    }

    if ion_name not in internal_data:
        raise ValueError(f"\"{ion_name}\" is invalid ion_name. Only NA, K or CL")

    template = "[ moleculetype ]\n"\
            "; name  nrexcl\n"\
            f"{internal_data[ion_name]['NAME']}  3\n\n"\
            f"[ atoms ]\n"\
            ";   nr   type  resnr residue  atom   cgnr     charge         mass\n"\
            f"   1     {internal_data[ion_name]['TYPE']}      1      {internal_data[ion_name]['RESIDUE']}    {internal_data[ion_name]['ATOM']}      1  {internal_data[ion_name]['CHARGE']}    {internal_data[ion_name]['MASS']}\n\n"
    return template

def add_ions_moleculetype(input_topology:PathLike, output_topology:PathLike):
    """Simple function to fix the topology file if the ions are not yet added
    as moleculetype

    Parameters
    ----------
    input_topology : PathLike
        Path to the input topology file
    output_topology : PathLike
        Path to the output topology file
    """

    molecules = get_molecule_names(input_topology, section='molecules')
    molecule_types = get_molecule_names(input_topology, section = 'moleculetype')
    ions = ['CL', 'NA', 'K']
    ions_moleculetype = ''
    for possible_ion in set(molecules) - set(molecule_types):
        if possible_ion in ions:
            ions_moleculetype += make_ion_moleculetype_section(possible_ion)
    
    with open(input_topology, "r") as topology_file:
        old_lines = topology_file.readlines()

    new_lines = []
    i = 0
    while i < len(old_lines):
        if '[ atomtypes ]' in old_lines[i]:
            new_lines.append(old_lines[i])
            i += 1
            while "[" not in old_lines[i]:
                new_lines.append(old_lines[i])
                i += 1
                if i >= len(old_lines): break
            # ffnonbonded.itp term
            new_lines.append(ions_moleculetype)
            new_lines.append(old_lines[i])
        else:
            new_lines.append(old_lines[i])
        i += 1
    with open(output_topology, 'w') as out:
        for line in new_lines:
            out.write(line)

def add_water_ions_param(input_topology:PathLike, output_topology:PathLike):
    """Add water and ion atom types to the main [ atomtypes ] section of the topology.
    Also [ moleculetype ], [ bonds ], [ angles ], [ settles ] and [ exclusions ]
    sections for the water molecule. This sections were taken from an example of BioSimSpace

    Parameters
    ----------
    input_topology : PathLike
        Path to the input topology file
    output_topology : PathLike
        Path to the output topology file
    """
    with open(input_topology, "r") as topology_file:
        old_lines = topology_file.readlines()

    new_lines = []
    i = 0
    while i < len(old_lines):
        if '[ atomtypes ]' in old_lines[i]:
            new_lines.append(old_lines[i])
            i += 1
            while "[" not in old_lines[i]:
                new_lines.append(old_lines[i])
                i += 1
                if i >= len(old_lines): break
            # ffnonbonded.itp term
            new_lines.append("OW           8      16.00    0.0000  A   3.15061e-01  6.36386e-01\n"\
                            "HW           1       1.008   0.0000  A   0.00000e+00  0.00000e+00\n"\
                            "; spc water - use only with spc.itp & settles\n"\
                            "OW_spc       8      15.9994  0.0000  A   3.16557e-01  6.50629e-01\n"\
                            "HW_spc       1       1.0080  0.0000  A   0.00000e+00  0.00000e+00\n"\
                            "K           19      39.10    0.0000  A   4.73602e-01  1.37235e-03\n"\
                            "Cl          17      35.45    0.0000  A   4.40104e-01  4.18400e-01\n"\
                            "Na          11      22.99    0.0000  A   3.32840e-01  1.15897e-02\n\n"
                            )
            new_lines.append(old_lines[i])

        elif '[ system ]' in old_lines[i]:
            new_lines.append("[ moleculetype ]\n"\
                "; name  nrexcl\n"\
                "SOL  3\n\n"\
                
                "[ atoms ]\n"\
                ";   nr   type  resnr residue  atom   cgnr     charge         mass\n"\
                "     1     OW      1     SOL    OW      1  -0.834000    16.000000\n"\
                "     2     HW      1     SOL   HW1      1   0.417000     1.008000\n"\
                "     3     HW      1     SOL   HW2      1   0.417000     1.008000\n\n"\
                
                "#ifdef FLEXIBLE\n"\
                "[ bonds ]\n"\
                ";   ai     aj  funct  parameters\n"\
                "     1      2      1  0.09572  50241\n"\
                "     1      3      1  0.09572  502416\n\n"\
                
                "[ angles ]\n"\
                ";   ai     aj     ak   funct   parameters\n"\
                "     2      1      3       1   104.52  628.02\n\n"\
                
                "#else\n\n"\

                "[ settles ]\n"\
                "; OW    funct   doh dhh\n"\
                "1       1       0.09572 0.15136\n\n"\
                
                "[ exclusions ]\n"\
                "1   2   3\n"\
                "2   1   3\n"\
                "3   1   2\n\n"\
                
                "#endif\n\n"\

                "[ system ]\n"
                )
        else:
            new_lines.append(old_lines[i])
        i += 1

    with open(output_topology, 'w') as out:
        for line in new_lines:
            out.write(line)

def fix_topology(input_topology: PathLike,
                 out_dir: PathLike,
                 exclusion_list:list = ["SOL", "NA", "CL", "MG", "ZN"],
                 f_xyz:tuple = (2500, 2500, 2500)):
    """It will go through input_topology, create the posres files for the identified molecules
    not in exclusion list, and finally add the corresponded include statements in the topology file.
    on out_topology_path you will have
    It will call, sequentially, to:

    #. :meth:`abfe.scripts.preparation.gmx_topology.get_molecule_names`
    #. :meth:`abfe.scripts.preparation.gmx_topology.make_posres_files`
    #. :meth:`abfe.scripts.preparation.gmx_topology.add_posres_section`

    In out_dir you will have:

    #. Fixed topology file name_of_input_topology_fix.top
    #. Posres itp files.

    Parameters
    ----------
    input_topology : PathLike
        The path of the input topology
    out_dir : PathLike
        Where the fixed and generated itp posres files will be written
    exclusion_list : list, optional
        Molecules names to do not take into account during the posres file generation, by default ["SOL", "NA", "CL", "MG", "ZN"]
    f_xyz : tuple
        The x, y, z components of the restraint force to be used. It could 
        be a float number of a string to be then defined on the mdp file, by default (2500, 2500, 2500)
    """
    name, ext = os.path.splitext(os.path.basename(input_topology))
    out_topology = os.path.join(out_dir, f"{name}_fix{ext}")
    molecules = list(set(get_molecule_names(input_topology)) - set(exclusion_list))
    make_posres_files(input_topology, out_dir = out_dir, molecules = molecules, f_xyz = f_xyz)
    add_posres_section(input_topology, molecules, out_file=out_topology)
    # In case that everything is fine, add_ions_moleculetype will not modify the topology
    add_ions_moleculetype(out_topology, out_topology)

def index_for_membrane_system(
        configuration_file:PathLike,
        ndxout:PathLike = "index.ndx",
        lignad_name:str = 'LIG',
        cofactor_name:str = None,
        cofactor_on_protein:bool = True):
    """Make the index file for membrane systems with SOLU, MEMB and SOLV. It uses gmx make_ndx and select internally.
    One examples selection that can be created with ligand_name = LIG; cofactor_name = COF and cofactor_on_protein = True is:
        #. "SOLU" group Protein or resname LIG or resname COF;
        #. "MEMB" ((group System and ! group Water_and_ions) and ! group Protein) and ! (resname LIG) and ! (resname COF);
        #. "SOLV" group Water_and_ions;


    Parameters
    ----------
    configuration_file : PathLike
        PDb or GRO file of the system.
    ndxout : PathLike
        Path to output the index file.
    lignad_name : str
        The residue name for the ligand in the configuration file, bt default LIG.
    cofactor_name : str
        The residue name for the cofactor in the configuration file, bt default None
    cofactor_on_protein : bool
        It only works if cofactor_name is provided. If True, the cofactor will be part of the protein and the lignad
        if False will be part of the solvent and ions, bt default True
    """
    tmpopt = tempfile.NamedTemporaryFile(suffix='.opt')
    tmpndx = tempfile.NamedTemporaryFile(suffix='.ndx')
    # Nice use of gmx select, see the use of the parenthesis
    sele_MEMB = f"\"MEMB\" ((group System and ! group Water_and_ions) and ! group Protein) and ! (resname {lignad_name})"
    sele_SOLU = f"\"SOLU\" group Protein or resname {lignad_name}"
    sele_SOLV = f"\"SOLV\" group Water_and_ions"
    if cofactor_name:
        sele_MEMB += f" and ! (resname {cofactor_name})"
        if cofactor_on_protein:
            sele_SOLU += f" or resname {cofactor_name}"
        else:
            sele_SOLV += f" or resname {cofactor_name}"          

    sele_SOLU += ";\n"
    sele_MEMB += ";\n"
    sele_SOLV += ";\n"
    print(sele_SOLU + sele_MEMB + sele_SOLV)
    with open(tmpopt.name, "w") as opt:
        opt.write(sele_SOLU + sele_MEMB + sele_SOLV)
    tools.run(f"""
                export GMX_MAXBACKUP=-1
                echo "q" | gmx make_ndx -f {configuration_file} -o {tmpndx.name}
                gmx select -s {configuration_file} -sf {tmpopt.name} -n {tmpndx.name} -on {ndxout}
                """)

    #deleting the line _f0_t0.000 in the file
    with open(ndxout, "r") as index:
        data = index.read()
        data = data.replace("_f0_t0.000","")
    with open(ndxout, "w") as index:
        index.write(data)

    tmpopt.close()
    tmpndx.close()


if __name__ == "__main__":...