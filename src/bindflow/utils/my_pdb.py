#!/usr/bin/env python

from string import ascii_uppercase
import numpy as np
import os, shutil

def backoff(file_path):
    """


    Parameters
    ----------
    file : TYPE string
        DESCRIPTION: The name or the path f     or the specific file

    Returns
    -------
    None.
    If the file already exist. it will made a back up to ./#{file}.{str(i)}#,
    Where i is an integer.
    """
    basname = os.path.basename(file_path)
    dirname = os.path.dirname(file_path)
    if os.path.exists(file_path):
        new_basname = basname
        i = 1
        while(os.path.exists(os.path.join(dirname, new_basname))):
            new_basname = f"./#{basname}.{str(i)}#"
            i += 1
        print(f"Back Off! I just backed up {file_path} to {os.path.join(dirname, new_basname)}")
        shutil.copy2(file_path, os.path.join(dirname, new_basname))


"""
Taken from https://charmm-gui.org/?doc=lecture&module=scientific&lesson=7
"""

class Atom:
    """
    https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
    """
    def __init__(self, line):
        self.code = "ATOM"
        self.serial = int(line[6:11])
        self.name = line[12:16].strip()
        self.altLoc = line[16]
        self.resName = line[17:21].strip()
        self.chainID = line[21]
        self.resSeq = int(line[22:26])
        self.iCode = line[26]
        self.x = float(line[30:38])
        self.y = float(line[38:46])
        self.z = float(line[46:54])
        self.occupancy = line[54:60].strip()
        self.tempFactor = line[60:66].strip()
        self.element = line[76:78].strip()
        self.charge = line[78:80].strip()
        if self.occupancy: self.occupancy = float(self.occupancy)
        if self.tempFactor: self.tempFactor = float(self.tempFactor)

    def __getitem__(self, key):
        return self.__dict__[key]


class Hetatm:
    """
    https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM
    """
    def __init__(self, line):
        self.code = "HETATM"
        self.serial = int(line[6:11])
        self.name = line[12:16].strip()
        self.altLoc = line[16]
        self.resName = line[17:21].strip()
        self.chainID = line[21]
        self.resSeq = int(line[22:26])
        self.iCode = line[26]
        self.x = float(line[30:38])
        self.y = float(line[38:46])
        self.z = float(line[46:54])
        self.occupancy = line[54:60].strip()
        self.tempFactor = line[60:66].strip()
        self.element = line[76:78].strip()
        self.charge = line[78:80].strip()
        if self.occupancy: self.occupancy = float(self.occupancy)
        if self.tempFactor: self.tempFactor = float(self.tempFactor)

    def __getitem__(self, key):
        return self.__dict__[key]


class CRYST1:
    """
    https://www.wwpdb.org/documentation/file-format-content/format33/sect8.html#CRYST1
    """
    def __init__(self, line):
        self.a = float(line[6:15])			    #Real(9.3)     a              a (Angstroms).
        self.b = float(line[15:24])			    #Real(9.3)     b              b (Angstroms).
        self.c = float(line[24:33])			    #Real(9.3)     c              c (Angstroms).
        self.alpha = float(line[33:40])			#Real(7.2)     alpha          alpha (degrees).
        self.beta = float(line[40:47])			#Real(7.2)     beta           beta (degrees).
        self.gamma = float(line[47:54])			#Real(7.2)     gamma          gamma (degrees).
        self.sGroup = line[55:66]			    #LString       sGroup         Space  group.
        try:
            self.z = int(line[66:70])			    #Integer       z              Z value.
        except:
               self.z = ""

    def __getitem__(self, key):
        return self.__dict__[key]
    def write(self):
        string = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%-12s%4s\n"%\
            (self.a,self.b,self.c,self.alpha,self.beta,self.gamma,self.sGroup,self.z)
        return string

class Residue:
    """
    This will convert a list of Atom objects that belongs and will create a
    residue class taken into account the information of the first atom.
    It suppose that all of them belong to the same residue
    """
    def __init__(self, atoms):
        self.resName = atoms[0].resName
        self.resSeq = atoms[0].resSeq
        self.chainID = atoms[0].chainID
        self.atoms = atoms

    def __getitem__(self, key):
        return self.__dict__[key]
    def __repr__(self) -> str:
        f"{self.__class__.__name__}(resName = {self.resName}, "\
        f"resSeq = {self.resSeq}, "\
        f"chainID = {self.chainID}, "
        f"NumAtoms = {len(self.atoms)})"

class PDB:
    def __init__(self, file:str):
        """The constructor of PDB

        Parameters
        ----------
        file : str
            The path of the pdb file
        """
        self.file = file
        self.title = ""
        self.cryst1 = CRYST1("0"*54+"  P 1 21 1     1")#Just to initialized to a generic cell, see the web referenced in CRYST1
        self.atoms = []
        self._read()

    def _read(self):
        """Reader of self.file
        """
        MODEL = None
        with open(self.file, 'r') as f:
            for line in f.readlines():
                if line.startswith('TITLE'): self.title = line.strip()
                elif line.startswith('CRYST1'): self.cryst1 = CRYST1(line)
                elif line.startswith('MODEL'): MODEL = int(line.split()[1])
                elif line.startswith('ATOM'):
                    atom = Atom(line)
                    atom.MODEL = MODEL
                    self.atoms.append(atom)
                elif line.startswith('HETATM'):
                    atom = Hetatm(line)
                    atom.MODEL = MODEL
                    self.atoms.append(atom)

    def get_chains(self) -> list[str]:
        """Get a list of chain identifiers

        Returns
        -------
        list[str]
            The list of chains
        """
        chains = set()
        for atom in self.atoms:
            chains.add(atom.chainID)
        return sorted(chains)


    def get_atoms(self, to_dict:bool = False):
        """Return a list of all atoms.

        Parameters
        ----------
        to_dict : bool, optional
            If to_dict is True, each atom is represented as a dictionary.
            Otherwise, a list of Atom objects is returned, by default True

        Returns
        -------
        list[Atoms] or list[dict]
            The list of atoms
        """
        if to_dict: return [x.__dict__ for x in self.atoms]
        else: return self.atoms


    def get_model(self, model_num, to_dict=True):
        """Return all atoms where MODEL == model_num"""
        model_atoms = [x for x in self.atoms if x.MODEL == model_num]
        if to_dict:
            return [atom.__dict__ for atom in model_atoms]
        else:
            return model_atoms


    def get_residues(self) -> list[Residue]:
        """Rearrange the data as residues function.

        Returns
        -------
        list[Residue]
            A list of Residues objects
        """
        residues = []
        i = 0
        while i < len(self.atoms):
            resAtoms = []
            resName = self.atoms[i].resName
            resSeq = self.atoms[i].resSeq
            chainID = self.atoms[i].chainID
            j = i
            while self.atoms[j].resName == resName and self.atoms[j].resSeq == resSeq and self.atoms[j].chainID == chainID:
                resAtoms.append(self.atoms[j])
                j += 1
                if j >= len(self.atoms): break
            residues.append(Residue(resAtoms))
            i = j

        return residues


    def edges(self):
        """
        print the edges of the structure
        """
        x_min = min(self.atoms, key=lambda item: item.x)
        y_min = min(self.atoms, key=lambda item: item.y)
        z_min = min(self.atoms, key=lambda item: item.z)

        x_max = max(self.atoms, key=lambda item: item.x)
        y_max = max(self.atoms, key=lambda item: item.y)
        z_max = max(self.atoms, key=lambda item: item.z)

        length = [round(x_max.x-x_min.x, 3), round(y_max.y-y_min.y, 3), round(z_max.z-z_min.z, 3)]

        print(f"""
        X axes:
            Length: {length[0]} Angstroms
            Bounding atoms: {x_max.name}_{x_max.serial} ({x_max.resName}_{x_max.resSeq}-{x_max.chainID}) ----- {x_min.name}_{x_min.serial} ({x_min.resName}_{x_min.resSeq}-{x_min.chainID})
        Y axes:
            Length: {length[1]} Angstroms
            Bounding atoms: {y_max.name}_{y_max.serial} ({y_max.resName}_{y_max.resSeq}-{y_max.chainID}) ----- {y_min.name}_{y_min.serial} ({y_min.resName}_{y_min.resSeq}-{y_min.chainID})
        Z axes:
            Length: {length[2]} Angstroms
            Bounding atoms: {z_max.name}_{z_max.serial} ({z_max.resName}_{z_max.resSeq}-{y_max.chainID}) ----- {z_min.name}_{z_min.serial} ({z_min.resName}_{z_min.resSeq}-{z_min.chainID})
        """)
        return tuple(length)

    def box_center(self, resName_selection:list[str] = None, H_atoms:bool = True):
        """Get the center of the box

        Parameters
        ----------
        resName_selection : list[str], optional
            A list of residues names, by default None
        H_atoms : bool, optional
            Add the hydrogens atoms to the calculation, by default True

        Returns
        -------
        tuple
            A tuple with the xyz coordinates of the box's center
        """
        if resName_selection:
            atoms = [atom for atom in self.atoms if atom.resName in resName_selection]
        else:
            atoms = self.atoms

        if not H_atoms:
            atoms = [atom for atom in atoms if atom.name[0] not in 'Hh']

        x_min = min(atoms, key=lambda item: item.x)
        y_min = min(atoms, key=lambda item: item.y)
        z_min = min(atoms, key=lambda item: item.z)

        x_max = max(atoms, key=lambda item: item.x)
        y_max = max(atoms, key=lambda item: item.y)
        z_max = max(atoms, key=lambda item: item.z)

        return round((x_max.x+x_min.x)/2,3), round((y_max.y+y_min.y)/2,3), round((z_max.z+z_min.z)/2,3)

    def centroid(self, resName_selection:list[str] = None, H_atoms:bool = True):
        """Get the centroid of the system

        Parameters
        ----------
        resName_selection : list[str], optional
            A list of residues names, by default None
        H_atoms : bool, optional
            Add the hydrogens atoms to the calculation, by default True

        Returns
        -------
        tuple
            A tuple with the xyz coordinates of the centroid
        """
        if resName_selection:
            atoms = [atom for atom in self.atoms if atom.resName in resName_selection]
        else:
            atoms = self.atoms

        if not H_atoms:
            atoms = [atom for atom in atoms if atom.name[0] not in 'Hh']

        centroid = (sum([np.array([atom.x,atom.y,atom.z]) for atom in atoms]) / len(atoms))
        centroid = [round(coord, 3) for coord in centroid]
        return tuple(centroid)


    def fix_chains(self, only_ATOM:bool = True):
        """Fix the chains on the PDB object.
        when there is only one amino acid per chain with the
        same residue number this will fail. It will assign chain A to all.

        Parameters
        ----------
        only_ATOM : bool, optional
            Use only the Atom, in False, HEATATM could be used, by default True
        """
        if only_ATOM:
            # TODO, check of this is right
            atoms = [atom for atom in self.atoms if atom.code == "ATOM"]
        else:
            atoms = self.atoms

        chain = 0
        for (i, atom) in enumerate(atoms):
            if i == 0:
                atom.chainID = ascii_uppercase[chain]
            elif atoms[i].resSeq >= atoms[i-1].resSeq:
                atom.chainID = ascii_uppercase[chain]
            else:
                chain += 1
            atom.chainID = ascii_uppercase[chain]

    def add_element(self):
        """
        This function add the columnd elemnt to the pdb data
        """
        for atom in self.atoms:
            if not atom.element:
                element = ""
                for char in atom.name:
                    if char.isalpha():
                        element = char.upper()
                        break
                atom.element = element

    def make_posres(self, heavy_atoms_only = True, posre_name = None, fcx_fcy_fcz = (1000,1000,1000)):
        """It create a position restraint itp file for GROMACS.
        The corresponded include statement is added to the topology.

        Parameters
        ----------
        heavy_atoms_only : bool, optional
            Only use the index of the heavy atoms, by default True
        posre_name : str, optional
            The name of the position restrain file, by default None
        fcx_fcy_fcz : tuple, optional
            The xyz component of the force constant, by default ('POSRES_LIG','POSRES_LIG','POSRES_LIG')
        """
        if not posre_name:
            posre_name = os.path.basename(self.file).split('.')[0]

        include_statement = "; Include Position restraint file for heavy atoms\n"\
                            "#ifdef POSRES\n"\
                            f"#include \"{posre_name}_posre.itp\"\n"\
                            "#endif"
        print(f"Add the following include statment to the GROMACS topology:\n {include_statement}")

        if heavy_atoms_only:
            atom_indexes = []
            for i, atom in enumerate(self.atoms):
                    if atom.name[0] != "h":
                        atom_indexes.append(i)
        else:
            atom_indexes = len(self.atoms)
        #Writing the posre file
        with open(f"{posre_name}_posre.itp", "w") as posre:
            posre.write("[ position_restraints ]\n")
            posre.write(f"{';i':<5}{'funct':^5}{'fcx':^20}{'fcy':^20}{'fcz':^20}\n")
            for atom_index in atom_indexes:
                posre.write(f"{atom_index:^5}{1:^5}{fcx_fcy_fcz[0]:^20}{fcx_fcy_fcz[1]:^20}{fcx_fcy_fcz[2]:^20}\n")

    def write(self, out_file = None, backup = True):
        """Write pdb file, all the modifications
        are taken into account.

        Parameters
        ----------
        outfile : str, optional
            The path to write the topology, by default None
        backup : bool, optional
            backup the original file in case of out_file has the same name or it was not provided, by default True
        """
        if out_file:
            to_work = out_file
        else:
            to_work = self.file
            if backup: backoff(to_work)
        with open(to_work,'w') as out:
            if self.title: out.write(self.title+'\n')
            if self.cryst1.a or self.cryst1.b or self.cryst1.c: #If the length of at least one vector is greater than 0
                out.write(self.cryst1.write())
            if self.atoms[0].MODEL != None: out.write("MODEL %6i\n"%(self.atoms[0].MODEL))


            # For the case that the protein have more than one chain, the TER
            # flag will be printed. If not (the protein has only one chain or
            # the chain is not declared) the TER flag will not printed

            #Because get_chains belongs to the class you need to call it like a self.function(), the method of the class
            chain_Ids = self.get_chains()
            if len(chain_Ids) > 1:
                last = self.atoms[0].chainID
                for (i,atom) in enumerate(self.atoms):

                    #INCREDIBLE, PropKa if the atoms name is shifted to the left dont understand!!!
                    #For that reason I need to change when the len is 3. NO IDEA!!!!
                    if len(atom.name) == 3:
                        string_name = f"{atom.name:>4}"
                    else:
                        string_name = f"{atom.name:^4}"

                    if atom.chainID != last:
                        out.write("TER\n")
                        last = atom.chainID

                    #This check is necesary for the case that some residues has a name greater than 3 letters (ASPH)
                    #This generate something like ASPHA, where A is the Chain identifier

                    out.write(f"{atom.code:<6}{atom.serial:5} {string_name}"
                              f"{atom.altLoc:1}{atom.resName:<4}{atom.chainID:1}"
                              f"{atom.resSeq:4}{atom.iCode:1}   {atom.x:8.3f}"
                              f"{atom.y:8.3f}{atom.z:8.3f}{atom.occupancy:6.2f}"
                              f"{atom.tempFactor:6.2f}           {atom.element:2}"
                              f"{atom.charge:2}\n")
                out.write("TER\n")#This is for the last TER flag
            else:
                for (i,atom) in enumerate(self.atoms):
                    #INCREDIBLE, PropKa if the atoms name is shifted to the left dont understand!!!
                    if len(atom.name) == 3:
                        string_name = f"{atom.name:>4}"
                    else:
                        string_name = f"{atom.name:^4}"

                    out.write(f"{atom.code:<6}{atom.serial:5} {string_name}"
                              f"{atom.altLoc:1}{atom.resName:<4}{atom.chainID:1}"
                              f"{atom.resSeq:4}{atom.iCode:1}   {atom.x:8.3f}"
                              f"{atom.y:8.3f}{atom.z:8.3f}{atom.occupancy:6.2f}"
                              f"{atom.tempFactor:6.2f}           {atom.element:2}"
                              f"{atom.charge:2}\n")


if __name__ == '__main__':
    #import sys
    #pdb = PDB(sys.argv[1])
    #pdb = PDB('PROTEIN.pdb')
    #print(pdb)
    #atoms = pdb.atoms
    #print(atoms[0].__dict__)
    #print(pdb.get_chains())
    #for atom in atoms: #Could be used to rectificate the pdb and put the chain information
    #    atom.chainID = "G"
    #print(pdb.get_chains())
    #pdb.cryst1.a = 500
    #print(pdb.cryst1.write())
    #for atom in atoms:
    #    if atom.resName == "LIA":
    #        atom.x = 0
    #        atom.y = 0
    #        atom.z = 0
    #pdb.write("j.pdb")
    #print(isinstance(pdb, PDB))
    #print(type(pdb).__name__)


    #print the name of each atom
    #for atom in atoms:
    #    print(atom.serial, atom.name, atom.x, atom.y, atom.z, atom.element)
    #lig_name = "HV6"
    #for atom in pdb.atoms:
    #    if atom.resName == lig_name:
    #        atom.resName = f"LI{atom.chainID}"
    #print(pdb.centroid(), pdb.box_center(), pdb.edges())
    #print(pdb.atoms[-1].__dict__)
    pass
