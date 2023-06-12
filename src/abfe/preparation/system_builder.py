#!/usr/bin/env python
import warnings
import glob
import os
import shutil
import argparse
import tempfile
import copy
from typing import Iterable
import tarfile
import logging
logger = logging.getLogger(__name__)

from abfe.home import home
from abfe.preparation.gmx_topology import fix_topology, add_water_ions_param, index_for_membrane_system
from abfe.utils.tools import run, PathLike

from toff import Parameterize

with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)
    import BioSimSpace as bss
    

from parmed.gromacs import GromacsTopologyFile, GromacsGroFile
from parmed.structure import Structure
from parmed.tools.actions import HMassRepartition
# from pdbfixer import PDBFixer
# from openmm.app import PDBFile

# TODO remove the BioSimSpace dependency and just use ParmED

def readParmEDMolecule(top_file:PathLike, gro_file:PathLike) -> Structure:
    """Read a gro and top GROMACS file and return
    a topology Structure

    Parameters
    ----------
    top_file : PathLike
        Path of the top file
    gro_file : PathLike
        Path of the gro file

    Returns
    -------
    Structure
        Structure with topologies, coordinates and box information
    """
    gmx_top = GromacsTopologyFile(top_file)
    gmx_gro = GromacsGroFile.parse(gro_file)

    # Add positions
    gmx_top.positions = gmx_gro.positions
    # Needed because .prmtop contains box info
    gmx_top.box = gmx_gro.box
    return gmx_top


def get_gmx_ff(ff_code:str, out_dir:PathLike = '.') -> PathLike:
    """Get GROMACS Force Field


    Parameters
    ----------
    ff_code : PathLike
        The identification of the gromacs force field.
        For now only: Slipids_2020 and amber99sb-star-ildn are supported.
    out_dir : PathLike, optional
        Where the file will be decompress, by default '.'
    """
    out_dir = os.path.abspath(out_dir)
    supported_ff = [
        'Slipids_2020',
        'amber99sb-star-ildn',
    ]
    if ff_code not in supported_ff:
        raise ValueError(f"ff_code = {ff_code} is not valid. Chose between: {supported_ff}")
    else:
        fname = os.path.join(home(dataDir='gmx_ff'), f'{ff_code}.ff.tar.gz')
    tar = tarfile.open(fname, "r:gz")
    tar.extractall(out_dir)
    tar.close()
    return os.path.join(out_dir,  f'{ff_code}.ff')


def system_combiner(**md_elements):
    """This function simply sum up all the elements provided 
    as keyword arguments.

    Returns
    -------
    object
        any Python object with the method sum implemented. In case elements
        that evaluate as False in Python will not be taken into account:
        E.g. False, 0, '', None

    Raises
    ------
    RuntimeError
        In case all the elements evaluate as False
    """
    if any(md_elements.values()):
        # md_system = sum(element for element in md_elements.values() if element) # it does not work with sum
        # Use copy to avoid inplace modifications
        for element in md_elements:
            if md_elements[element]:
                try:
                    md_system += copy.copy(md_elements[element])
                except NameError:
                    md_system = copy.copy(md_elements[element])
    else:
        raise RuntimeError(f"\t* system_combiner failed with the inputs: {md_elements}")
    print(f"\t* The system was constructed as fallow: {' + '.join([key for key in md_elements if md_elements[key]])}")
    return md_system

def parmed_solvate(
        parmed_system:Structure,
        bt:str = "triclinic",
        box:list[float] = None, # in nm, PDB are in Angstroms
        angles:list[float] = None, # In degrees (90,90,60) for membrane systems
        d:float = None,
        c:bool = False,
        cs:str = "spc216",
        pname:str = "K",
        nname:str = "CL",
        rmin:float = 1.0,
        out_dir:PathLike = '.'):
    
    # TODO, ion concentration is not used, only neutralize
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    cwd = os.getcwd()
    with tempfile.TemporaryDirectory(prefix='solv_', dir=cwd) as tmp:
        os.chdir(tmp)
        parmed_system.save('non_water_ions_system.top', overwrite = True)
        parmed_system.save('system.gro', overwrite = True)


        editconf_cmd = f"gmx editconf -f system.gro -o system_box_corrected.gro -bt {bt}"
        if box: editconf_cmd += f" -box {' '.join([str(i) for i in box])}"
        if angles: editconf_cmd += f" -angles {' '.join([str(i) for i in angles])}"
        if d: editconf_cmd += f" -d {d}"
        if c: editconf_cmd += " -c"
        
        # First write an mdp file.
        with open("ions.mdp", "w") as file:
            file.write("; Neighbor searching\n"\
            "cutoff-scheme           = Verlet\n"\
            "rlist                   = 1.1\n"\
            "pbc                     = xyz\n"\
            "verlet-buffer-tolerance = -1\n"\
            "\n; Electrostatics\n"\
            "coulombtype             = cut-off\n"\
            "\n; VdW\n"\
            "rvdw                    = 1.0\n")
        

        run(f"""
            export GMX_MAXBACKUP=-1
            {editconf_cmd}
            gmx solvate -cp system_box_corrected.gro -p non_water_ions_system.top -cs {cs} -o system_solvated.gro
        """)
        add_water_ions_param('non_water_ions_system.top', 'system.top')
        run(f"""    
            gmx grompp -f ions.mdp -c system_solvated.gro -p system.top -o ions.tpr
            echo "SOL" | gmx genion -s ions.tpr -p system.top -o system_solvated_ions.gro -neutral -pname {pname} -nname {nname} -rmin {rmin}
        """)
        final_top = os.path.abspath('system.top')
        final_gro = os.path.abspath('system_solvated_ions.gro')
        os.chdir(cwd)
        shutil.copy(final_top, os.path.join(out_dir, 'solvated.top'))
        shutil.copy(final_gro, os.path.join(out_dir, 'solvated.gro'))




# TODO, check what is the type of the bss_systems to add it as a HintType
def bss_solvate(bss_system:object, out_dir:PathLike = '.', vectors:Iterable[float] = None, angles:Iterable[float] = None):
    """Solvate and add ions to the system, if vectors and angles are not provided,
    the system will be solvated as a truncated octahedron with a padding of 15 Angstroms.

    Parameters
    ----------
    bss_system : object
        The BSS system to solvate
    out_dir : PathLike, optional
        Where the files will be written: solvated.gro, solvated.top, by default '.'
    vectors : Iterable[float], optional
        This is the vectors of the bos in ANGSTROMS. It is important that the provided vector has the correct units, by default None
    angles : Iterable[float], optional
        This is the angles between the components of the vector in DEGREES. It is important that the provided vector has the correct units, by default None

    Raises
    ------
    ValueError
        if vectors does not have three elements
    ValueError
        if angles does not have three elements
    """
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if isinstance(bss_system, Structure):
        pass
    else:
        if vectors and angles:
            if len(vectors) != 3:
                raise ValueError(f"'vectors' must be a iterable of three values: [a,b,c] in Angstrom. Provided: vectors = {vectors}")
            elif len(angles) != 3:
                raise ValueError(f"'angles' must be a iterable of three values: [alpha,beta,gamma] in Degree. Provided: angles = {angles}")
            vectors = [bss.Types.Length(proj, 'angstrom') for proj in vectors]
            angles = [bss.Types.Angle(angle, 'degree') for angle in angles]
        else:
            box_min, box_max = bss_system.getAxisAlignedBoundingBox()
            box_size = [y - x for x, y in zip(box_min, box_max)]
            padding = 15 * bss.Units.Length.angstrom
            box_length = (max(box_size) + 1.5 * padding)
            vectors, angles = bss.Box.truncatedOctahedron(box_length.value() * bss.Units.Length.angstrom)
        
        # TODO, ion concentration is not used, only neutralize
        solvated = bss.Solvent.tip3p(bss_system, box=vectors, angles=angles)#, ion_conc=)
        
        cwd = os.getcwd()
        os.chdir(out_dir)
        bss.IO.saveMolecules('solvated', solvated, ["GroTop", "Gro87"])
        os.chdir(cwd)

def make_abfe_dir(out_dir:PathLike, ligand_dir:PathLike, sys_dir:PathLike):
    """A copy and paste function to create the structure of the abfe directory

    Parameters
    ----------
    out_dir : PathLike
        Where the complex and the ligand systems will be created
    ligand_dir : PathLike
        Origin of the ligand inputs configuration and topologies files
    sys_dir : PathLike
        Origin of the complex inputs configuration and topologies files
    """
    complex_out = os.path.join(out_dir, "complex")
    ligand_out = os.path.join(out_dir, "ligand")
    if (not os.path.exists(complex_out)): os.makedirs(complex_out)
    if (not os.path.exists(ligand_out)): os.makedirs(ligand_out)

    for itp_ndx_file in glob.glob(os.path.join(ligand_dir, "*.itp")) + glob.glob(os.path.join(ligand_dir, "*.ndx")):
        shutil.copy(src=itp_ndx_file, dst=ligand_out)

    shutil.copyfile(src=os.path.join(ligand_dir, "solvated.gro"), dst=os.path.join(ligand_out, "ligand.gro"))
    shutil.copyfile(src=os.path.join(ligand_dir, "solvated_fix.top"), dst=os.path.join(ligand_out, "ligand.top"))

    for itp_ndx_file in glob.glob(os.path.join(sys_dir, "*.itp")) + glob.glob(os.path.join(sys_dir, "*.ndx")):
        shutil.copy(src=itp_ndx_file, dst=complex_out)

    shutil.copyfile(src=os.path.join(sys_dir, "solvated.gro"), dst=os.path.join(complex_out, "complex.gro"))
    # The last one in be copy, this will be used in the snake rule
    shutil.copyfile(src=os.path.join(sys_dir, "solvated_fix.top"), dst=os.path.join(complex_out, "complex.top"))


class CRYST1:
    """
    https://www.wwpdb.org/documentation/file-format-content/format33/sect8.html#CRYST1
    """
    def __init__(self, line = None):
        if line:
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
            self.__is_init = True
        else:
            self.__is_init = False
    
    def from_pdb(self, file:PathLike):
        """Initialize the class from a pdb file

        Parameters
        ----------
        file : PathLike
            The PDB file
        """
        with open(file, 'r') as f:
            for line in f.readlines():
                if line.startswith('CRYST1'):
                    self.__init__(line)
                    self.__is_init = True
                    break
        if not self.__is_init:
            warnings.warn('from_pdb was not able to initialize {self.__class__.__name__}')
    
    def get_bss_vectors(self) -> tuple[bss.Types.Length]:
        """get BioSimSpace vectors from the CRYST1 information

        Returns
        -------
        tuple[bss.Types.Length]
            BioSimSpace box (a, b, c)
        """
        return bss.Types.Length(self.a, 'angstrom'), bss.Types.Length(self.b, 'angstrom'), bss.Types.Length(self.c, 'angstrom')

    def get_bss_angles(self) -> tuple[bss.Types.Angle]:
        """get BioSimSpace angles from the CRYST1 information

        Returns
        -------
        tuple[bss.Types.Angle]
            BioSimSpace box (alpha, beta, gamma)
        """
        return bss.Types.Angle(self.alpha, 'degree'), bss.Types.Angle(self.beta, 'degree'), bss.Types.Angle(self.gamma, 'degree')

    def __getitem__(self, key):
        return self.__dict__[key]
    
    def string(self):
        string_repr = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%-12s%4s\n"%\
            (self.a,self.b,self.c,self.alpha,self.beta,self.gamma,self.sGroup,self.z)
        return string_repr
    
    def __repr__(self):
        return self.string()

class MakeInputs:
    def __init__(self,
            protein_pdb:PathLike = None,
            membrane_pdb:PathLike = None,
            cofactor_mol:PathLike = None,
            cofactor_on_protein:bool = True,
            hmr_factor:float = None,
            builder_dir:PathLike = 'builder'):
        """This class is used for building the systems for ABFE calculation.
        It will create the necessary topology and configuration files, as well the
        correct directory trees.

        Parameters
        ----------
        protein_pdb : PathLike, optional
            Path of the protein PDB file. Topology will be generated from amber99sb-ildn, by default None
        membrane_pdb : PathLike, optional
            Path of the membrane PDB file.
            If provided, the PDB must have a correct definition of the CRYST1,
            that information will be used for the solvation step,
            Topology will be generated from SLipid2020. by default None
        cofactor_mol : PathLike, optional
            Path of the mol cofactor file. Topology will be generated from OpenFF, by default None
        cofactor_on_protein : bool
            It is used during the index generation for membrane systems. It only works if cofactor_mol is provided.
            If True, the cofactor will be part of the protein and the ligand
            if False will be part of the solvent and ions, bt default True
        hmr_factor : float, optional
            The Hydrogen Mass Factor to use, by default None
        builder_dir : PathLike, optional
            Where all the building files. After completion you can safely remove calling the method clean, by default builder
        """
        self.protein_pdb = protein_pdb
        self.membrane_pdb = membrane_pdb
        self.cofactor_mol = cofactor_mol
        self.cofactor_on_protein = cofactor_on_protein
        self.hmr_factor = hmr_factor
        self.wd = os.path.abspath(builder_dir)
        self.__self_was_called = False

        # Initialize vectors and angles based on the information of the PDB only if a membrane system
        if self.membrane_pdb:
            cryst_info = CRYST1()
            cryst_info.from_pdb(self.membrane_pdb)
            self.vectors = (cryst_info.a/10, cryst_info.b/10, cryst_info.c/10) # Must convert form Angstrom to nm
            self.angles = (cryst_info.alpha, cryst_info.beta, cryst_info.gamma)
            
            logger.info(f"This is a membrane system. Crystal information was taken from {self.membrane_pdb} and it will be used for solvating the system as a GROMACS triclinic box: \n\t\t{cryst_info}")
        else:
            self.vectors, self.angles  = None, None


    def openff_process(self, mol_file:PathLike, name:str="MOL", safe_naming_prefix:str = None):
        """Get parameters for small molecules: ligands, cofactors, ...

        Parameters
        ----------
        mol_file : PathLike
            The path where the molecule is
        name : str, optional
            Name to give, by default "MOL"
        safe_naming_prefix : str, optional
            This is used to be sure that there will not happen any naming conflict in hte topologies, by default None

        Returns
        -------
        object
            The BioSimSpace system
        """

        if mol_file:
            print(f'\t\t- Getting OpenFF parameters for: {mol_file}')
        else:
            return None
        
        parameterizer = Parameterize(
            force_field_code = 'openff_unconstrained-2.0.0.offxml',
            ext_types = ['top', 'gro'],
            hmr_factor = self.hmr_factor,
            overwrite = True,
            safe_naming_prefix = safe_naming_prefix,
            out_dir = self.wd,
        )
        # Actually you can pass to parameterize Chem.rdchem.Mol, *.inchi, *.smi, *.mol, *.mol2, *.sdf
        parameterizer(input_mol = mol_file,mol_resi_name = name)

        top_file = os.path.join(self.wd, f"{name}.top")
        gro_file = os.path.join(self.wd, f"{name}.gro")
        
        if self.membrane_pdb:
            parmed_system = readParmEDMolecule(top_file=top_file, gro_file = gro_file)
            return parmed_system
        else:
            bss_system = bss.IO.readMolecules([top_file, gro_file])
            return bss_system

    def gmx_process(self, pdb_file:PathLike, pH:float = 7.0, is_membrane = False):
        """Used to process those biomolecules compatibles with amber99sb-ildn (protein, DNA, ..)
        and membrane compatibles with Slipid2020

        Parameters
        ----------
        pdb_file : PathLike
            Path to the file
        pH : float, optional
            pH for protonation (not working at the moment), by default 7.0
        is_membrane : bool, optional
            If True, Slipid2020 will be used instead of amber99sb-ildn, by default False

        Returns
        -------
        object
            The BioSimSpace system
        """

        if pdb_file:
            if is_membrane:
                print(f'\t\t- Getting Slipid2020 parameters for: {pdb_file}')
            else:
                print(f'\t\t- Getting amber99sb-ildn parameters for: {pdb_file}')
        else:
            return None
        
        name, _ = os.path.splitext(os.path.basename(pdb_file))

        # TODO, chake what is going wrong, and use this kind of code, much better than call from the command line.
        # fixer = PDBFixer(filename=pdb_file)
        # if not is_membrane:
        #     fixer.findMissingResidues()
        #     fixer.findNonstandardResidues()
        #     fixer.replaceNonstandardResidues()
        #     fixer.removeHeterogens(True)
        # fixer.findMissingAtoms()
        # if not is_membrane:
        #     fixer.addMissingAtoms()
        # fixer.addMissingHydrogens(pH)
        # with open(fixed_pdb, 'w') as f:
        #     PDBFile.writeFile(fixer.topology, fixer.positions, f)

        env_prefix = os.environ["CONDA_PREFIX"]
        # run(f"{env_prefix}/bin/pdbfixer {pdb_file} --output={fixed_pdb} --add-atoms=all --replace-nonstandard")

        gro_out = os.path.join(self.wd, f'{name}.gro')
        top_out = os.path.join(self.wd, f'{name}.top')
        posre_out = os.path.join(self.wd, f'{name}_posre.itp')

        if is_membrane:
            pdb_file = os.path.abspath(pdb_file)
            cwd = os.getcwd()
            get_gmx_ff('Slipids_2020', out_dir=self.wd)
            os.chdir(self.wd)
            run(f"gmx pdb2gmx -f {pdb_file} -ff Slipids_2020 -water none -o {gro_out} -p {top_out} -i {posre_out}")
            os.chdir(cwd)
        else:
            fixed_pdb = os.path.join(self.wd,f"{name}_fixed.pdb")
            run(f"{env_prefix}/bin/pdbfixer {pdb_file} --output={fixed_pdb} --add-atoms=all --replace-nonstandard")
            run(f"gmx pdb2gmx -f {fixed_pdb} -merge all -ff amber99sb-ildn -water tip3p -o {gro_out} -p {top_out} -i {posre_out} -ignh")

        if self.membrane_pdb:
            # BioSimSpace does not know how to deal with SLipid2020
            system = readParmEDMolecule(top_file=top_out, gro_file=gro_out)
            if self.hmr_factor:
                HMassRepartition(system, self.hmr_factor).execute()
            system.write(os.path.join(self.wd, f'{name}_final.top'))
        else:
            system = bss.IO.readMolecules([gro_out,top_out])
            if self.hmr_factor:
                system.repartitionHydrogenMass(factor=self.hmr_factor, water="no")
            cwd = os.getcwd()
            os.chdir(self.wd)
            bss.IO.saveMolecules(f"{name}_final", system, ["GroTop"])
            system = bss.IO.readMolecules([f"{name}_final.top", f'{name}.gro'])
            os.chdir(cwd)
        
        return system

    def make_system(self, ligand_mol:PathLike):
        """Create self.sys_ligand, self.sys_cofactor, self.sys_protein, self.sys_membrane
        and self.md_system (the combination of the available components). In case
        that the class was already called, it will be assumed that self.sys_cofactor, self.sys_protein, self.sys_membrane
        ere already calculated, only self.sys_ligand will be updated as well self.md_system

        Parameters
        ----------
        ligand_mol : PathLike
            Path of the ligand mol file
        """
        print("\t* Processing system components")
        self.sys_ligand = self.openff_process(
            mol_file = ligand_mol,
            name="LIG",
            safe_naming_prefix='x')
        
        # Only if the class has not yet called the full build will be carry out.
        if self.__self_was_called:
            print(f"\t\t- Reusing components from cache")
        else:
            self.sys_cofactor = self.openff_process(
                mol_file = self.cofactor_mol,
                name="COF",
                safe_naming_prefix='z')
            self.sys_protein = self.gmx_process(pdb_file = self.protein_pdb)
            self.sys_membrane = self.gmx_process(pdb_file = self.membrane_pdb, is_membrane = True)
        print("\t\t- Merging Components")
        # Cofactor at the end in case is a water molecule, not complains from GROMACS after solvation
        self.md_system = system_combiner(protein=self.sys_protein, membrane=self.sys_membrane, ligand=self.sys_ligand, cofactor=self.sys_cofactor)

    def clean(self):
        """Small cleaner, the intermediate steps saved on builder_dir will be deleted
        """
        try:
            shutil.rmtree(self.wd)
        except FileNotFoundError:
            pass

    def __call__(self, ligand_mol:PathLike, out_dir = 'abfe'):
        """The call implementation. It identify if it is needed to build
        all the components of the systems,
        In case that the class was already called, it will assume that all the components of the system,
        with the exception of the ligand, were already builded. This is useful to call the class
        on several ligands that share the same components: protein, membrane and cofactor

        Parameters
        ----------
        ligand_mol : PathLike
            Path of the mol ligand file. Topology will be generated from OpenFF
        out_dir : str, optional
            Where you would like to export the generated files, by default 'abfe'
        """
        print('\n\t'+42*"-")
        print(f"Processing ligand: {ligand_mol}")
        # Update (on multiple calls) or just create the out_dir (first call)
        self.out_dir = out_dir
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        
        # Construct MD system:
        self.make_system(ligand_mol)
        system_dir = os.path.join(self.wd, 'system')
        ligand_dir = os.path.join(self.wd, 'ligand')

        print("\t* Solvating:")
        if self.membrane_pdb:
            print("\t\t- Complex in: ", system_dir)
            parmed_solvate(self.md_system, bt='triclinic', box = self.vectors,angles=self.angles, out_dir=system_dir)
            print("\t\t- Ligand in: ", system_dir)
            parmed_solvate(self.sys_ligand, bt='octahedron', d = 1.5, out_dir=ligand_dir)
        else:
            print("\t\t- Complex in: ", system_dir)
            bss_solvate(self.md_system, out_dir=system_dir)
            print("\t\t- Ligand in: ", system_dir)
            bss_solvate(self.sys_ligand, out_dir=ligand_dir)

        print("\t* Fixing topologies")
        # TODO, I am here!! in case of membrane system, create index file and specific restraints files
        
        # Set the propers constraints depending on the system
        if self.membrane_pdb:
            # TODO: this is the easiest way to implement the position restraints changing the restraints
            # during different steps. However, we are not using different restraints for different molecules
            # what might be needed for some systems. That will take some coding in order to identify the molecules
            f_xyz = 3*['POSRES_DYNAMIC']
        else:
            f_xyz = 3*[2500]
        fix_topology(input_topology=os.path.join(system_dir,'solvated.top'), out_dir=system_dir, f_xyz=f_xyz)
        # For the ligand, I will keep the original f_xyz
        fix_topology(input_topology=os.path.join(ligand_dir,'solvated.top'), out_dir=ligand_dir, f_xyz=3*[2500])
        
        # Make index file in case of membrane systems
        if self.membrane_pdb:
            index_for_membrane_system(
                configuration_file = os.path.join(system_dir, "solvated.gro"),
                ndxout = os.path.join(system_dir, "index.ndx"),
                lignad_name = 'LIG',
                cofactor_name = 'COF' if self.cofactor_mol else None,
                cofactor_on_protein=self.cofactor_on_protein,
            )
        
        # Construct ABFE system:
        print(f"\t* Final build of ABFE directory on: {self.out_dir}")
        make_abfe_dir(out_dir=self.out_dir, ligand_dir=ligand_dir, sys_dir=system_dir)

        # Change state
        self.__self_was_called = True
        print("\n\t--------- Calculation Completed ----------\n")

def __system_builder_cmd():
    """
    Command line implementation for :meth:`abfe.system_builder.MakeInputs`
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        help='The directory where all the ligans are in the format of .mol',
        dest = 'ligand_mols_path',
        type = str,
    )
    parser.add_argument(
        '--protein_pdb',
        help = 'Protein pdb file, by default %(default)s',
        dest = 'protein_pdb',
        default = None,
        type = str,
    )
    parser.add_argument(
        '--membrane_pdb',
        help='Membrane pdb file, by default %(default)s',
        dest = 'membrane_pdb',
        default = None,
        type = str,
    )
    parser.add_argument(
        '--cofactor_mol',
        help='Cofactor mol file, by default %(default)s',
        dest = 'cofactor_mol',
        default = None,
        type = str,
    )
    parser.add_argument(
        '--hmr_factor',
        help='Hydrogen Mass Repartition factor, by default %(default)s',
        dest = 'hmr_factor',
        default = 3.0,
        type = float,
    )
    parser.add_argument(
        '--builder_dir',
        help='Where all the building files. After completion you can safely remove calling the method clean, by default %(default)s',
        dest = 'builder_dir',
        default = 'abfe',
        type = str,
    )
    parser.add_argument(
        '--out_dir',
        help='The directory where the build systems will be output, by default %(default)s',
        dest = 'out_dir',
        default = 'abfe',
        type = str,
    )
    args = parser.parse_args()
    builder = MakeInputs(
        protein_pdb=args.protein_pdb,
        membrane_pdb=args.membrane_pdb,
        cofactor_mol=args.cofactor_mol,
        hmr_factor=args.hmr_factor,
        builder_dir = args.builder_dir,
    )
    for ligand in glob.glob(os.path.join(args.ligand_mols_path, '*.mol')):
        name, _ = os.path.splitext(os.path.basename(ligand))
        builder(
            ligand_mol=ligand,
            out_dir= os.path.join(args.out_dir, name))
    builder.clean()


#############################################################################################

if __name__ == "__main__": ...