#!/user/bin python
import os
from typing import Union, Iterable
PathLike = Union[os.PathLike, str, bytes]
from dataclasses import dataclass

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
        The list of name of the molecules for which the topology section will be add
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

def make_posres_files(input_topology:PathLike, molecules:Iterable[str], out_dir:PathLike):
    """Make a position restraint file out of input_topology for all the molecules specified
    on molecules. The force constant used will be: fx = 2500; fy = 2500; fz = 2500
    for all atoms (included hydrogens, TODO: Check if this is what is expected, or only on heavy atoms)

    Parameters
    ----------
    input_topology : PathLike
        The path of the input topology
    molecules : Iterable[str]
        The list of name of the molecules for which the posres file will be created
    out_dir : PathLike
        The path where the posres files will be written
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
                                if float(top_lines[j].split()[7]) > 3:
                                    posres_str = f"{top_lines[j].split()[0]} 1 2500 2500 2500\n"
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

def fix_topology(input_topology: PathLike, out_dir: PathLike, exclusion_list:list = ["SOL", "NA", "CL", "MG", "ZN"]):
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
    """
    name, ext = os.path.splitext(os.path.basename(input_topology))
    out_topology = os.path.join(out_dir, f"{name}_fix{ext}")
    molecules = list(set(get_molecule_names(input_topology)) - set(exclusion_list))
    make_posres_files(input_topology, out_dir = out_dir, molecules = molecules)
    add_posres_section(input_topology, molecules, out_file=out_topology)
    # In case that everything is fine, add_ions_moleculetype will not modify the topology
    add_ions_moleculetype(out_topology, out_topology)


if __name__ == "__main__":...
    # input_topology = '/home/users/alejandro/GIT/ABFE_workflow/examples/prepearing_system/abfe/abc/input/complex/complex.top'
    
    # add_ions_moleculetype(input_topology, 'top.top')


    # # input_topology_path = sys.argv[1]
    # # out_topology_path = sys.argv[2]
    # exclusion_list = ["SOL", "NA", "CL", "MG", "ZN"]
    # fix_topology(
    #     input_topology = '/home/users/alejandro/GIT/ABFE_workflow/examples/prepearing_system/abfe/abc/input/complex/complex.top',
    #     out_dir='test',
    #     exclusion_list=exclusion_list,

    # )

    # # fix_topology(input_topology_path=input_topology_path, out_topology_path=out_topology_path, exclusion_list=exclusion_list)