#!/usr/bin/env python
import copy
import glob
import logging
import os
import shutil
import tarfile
import warnings
from typing import Union

from parmed.gromacs import GromacsGroFile, GromacsTopologyFile
from parmed.structure import Structure
from parmed.tools.actions import HMassRepartition
from toff import Parameterize

from abfe.home import home
from abfe.preparation import solvent
from abfe.utils.tools import PathLike, recursive_update_dict, run

# from pdbfixer import PDBFixer
# from
logger = logging.getLogger(__name__)


def readParmEDMolecule(top_file: PathLike, gro_file: PathLike) -> Structure:
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


def get_gmx_ff(ff_code: str, out_dir: PathLike = '.') -> PathLike:
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
    logger.info(f" The system was constructed as follows: {' + '.join([key for key in md_elements if md_elements[key]])}")
    return md_system


def make_abfe_dir(out_dir: PathLike, ligand_dir: PathLike, sys_dir: PathLike):
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
    if (not os.path.exists(complex_out)):
        os.makedirs(complex_out)
    if (not os.path.exists(ligand_out)):
        os.makedirs(ligand_out)

    for itp_ndx_file in glob.glob(os.path.join(ligand_dir, "*.itp")) + glob.glob(os.path.join(ligand_dir, "*.ndx")):
        shutil.copy(src=itp_ndx_file, dst=ligand_out)

    shutil.copyfile(src=os.path.join(ligand_dir, "solvated.gro"), dst=os.path.join(ligand_out, "ligand.gro"))
    shutil.copyfile(src=os.path.join(ligand_dir, "solvated.top"), dst=os.path.join(ligand_out, "ligand.top"))

    for itp_ndx_file in glob.glob(os.path.join(sys_dir, "*.itp")) + glob.glob(os.path.join(sys_dir, "*.ndx")):
        shutil.copy(src=itp_ndx_file, dst=complex_out)

    shutil.copyfile(src=os.path.join(sys_dir, "solvated.gro"), dst=os.path.join(complex_out, "complex.gro"))
    # The last one in be copy, this will be used in the snake rule
    shutil.copyfile(src=os.path.join(sys_dir, "solvated.top"), dst=os.path.join(complex_out, "complex.top"))


class CRYST1:
    """https://www.wwpdb.org/documentation/file-format-content/format33/sect8.html#CRYST1
    """
    def __init__(self, line: str = None):
        """The constructor

        Parameters
        ----------
        line : str, optional
            a line to process, by default None
        """
        if line:
            self.a = float(line[6:15])			    # Real(9.3)     a              a (Angstroms).
            self.b = float(line[15:24])			    # Real(9.3)     b              b (Angstroms).
            self.c = float(line[24:33])			    # Real(9.3)     c              c (Angstroms).
            self.alpha = float(line[33:40])			# Real(7.2)     alpha          alpha (degrees).
            self.beta = float(line[40:47])			# Real(7.2)     beta           beta (degrees).
            self.gamma = float(line[47:54])			# Real(7.2)     gamma          gamma (degrees).
            self.sGroup = line[55:66]			    # LString       sGroup         Space  group.
            try:
                self.z = int(line[66:70])			# Integer       z              Z value.
            except:
                self.z = ""
            self.__is_init = True
        else:
            self.__is_init = False

    def from_pdb(self, file: PathLike):
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

    def __getitem__(self, key):
        return self.__dict__[key]

    def string(self):
        string_repr = "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f%-12s%4s\n" %\
            (self.a, self.b, self.c, self.alpha, self.beta, self.gamma, self.sGroup, self.z)
        return string_repr

    def __repr__(self):
        return self.string()


class MakeInputs:
    """This class is used for building the systems for ABFE calculation.
    It will create the necessary topology and configuration files, as well the
    correct directory trees.
    """
    def __init__(self, protein: dict = None, membrane: dict = None, cofactor: dict = None,
                 cofactor_on_protein: bool = True, water_model: str = 'amber/tip3p',
                 custom_ff_path: Union[None, PathLike] = None, hmr_factor: Union[float, None] = None,
                 builder_dir: PathLike = 'builder'):
        """Constructor

        Parameters
        ----------
        protein : PathLike, optional
            This is a dictionary with the following information for the protein:

                * conf -> The path of the protein PDB/GRO file [mandatory]

                * top -> GROMACS topology [optional], by default None.
                Should be a single file topology with all the force field
                information and without the position restraint included. However, in case,
                you need to use an include statement such as:

                    include "./charmm36-jul2022.ff/forcefield.itp"

                You must change the statement to:

                    include "charmm36-jul2022.ff/forcefield.itp"

                And copy the charmm36-jul2022.ff to custom_ff_path and set this parameter accordingly. If not
                you may get some errors about files not founded. The force field directory
                must end with the suffix ".ff".

                * ff
                    * code -> GMX force field code [optional], by default amber99sb-ildn
                    You can use your custom force field, but custom_ff_path must be provided

        membrane : PathLike, optional
            This is a dictionary with the following information for the membrane:

                * conf -> The path of the membrane PDB file [mandatory]. If provided, the PDB must have a
                correct definition of the CRYST1. This information will be used for the solvation step

                * top -> GROMACS topology [optional], by default None.
                Should be a single file topology with all the force field
                information and without the position restraint included. However, in case,
                you need to use an include statement such as:

                    include "./amber-lipids14.ff/forcefield.itp"

                You must change the statement to:

                    include "/amber-lipids14.ff/forcefield.itp"

                And copy theamber-lipids14.ff to custom_ff_path and set this parameter accordingly. If not
                You may get some errors about files not founded. The force field directory
                must end with the suffix ".ff".

                * ff

                    * code -> GMX force field code [optional], by default Slipids_2020
                    You can use yoru custom force field, but custom_ff_path must be provided

        cofactor : dict, optional

            This is a dictionary with the following information for the cofactor:

                * conf -> The path of the small molecule file [mandatory]

                * top -> GROMACS topology [optional]. Must be a single file topology with all the
                force field information and without the position restraint included, by default None

                * ff

                    * code -> OpenFF code [optional], by default openff_unconstrained-2.0.0.offxml

        cofactor_on_protein : bool
            It is used during the index generation for membrane systems. It only works if cofactor_mol is provided.
            If True, the cofactor will be part of the protein and the ligand
            if False will be part of the solvent and ions, by default True
        hmr_factor : float, optional
            The Hydrogen Mass Factor to use, by default None
        water_model : str, optional
            The water force field to use, by default amber/tip3p
        custom_ff_path Union[None, PathLike], optional
            All the custom force field must be in this directory. The class will set:

                os.environ["GMXLIB"] = os.path.abspath(custom_ff_path)

        builder_dir : PathLike, optional
            Where all the building files. After completion you can safely remove calling the method clean, by default builder
        """
        self.protein = protein
        self.membrane = membrane
        self.cofactor = cofactor
        self.cofactor_on_protein = cofactor_on_protein
        self.hmr_factor = hmr_factor
        self.water_model = water_model
        self.wd = os.path.abspath(builder_dir)
        os.makedirs(self.wd, exist_ok=True)
        self.__self_was_called = False

        # Setting environmental variable for user custom force field:
        if custom_ff_path:
            self.custom_ff_path = os.path.abspath(custom_ff_path)
            os.environ["GMXLIB"] = self.custom_ff_path
        else:
            self.custom_ff_path = None

        # Initialize vectors and angles based on the information of the PDB only if a membrane system
        if self.membrane:
            cryst_info = CRYST1()
            cryst_info.from_pdb(self.membrane['conf'])
            self.vectors = (cryst_info.a/10, cryst_info.b/10, cryst_info.c/10)  # Must convert form Angstrom to nm
            self.angles = (cryst_info.alpha, cryst_info.beta, cryst_info.gamma)

            logger.info(f"This is a membrane system. Crystal information was taken from the configuration of: "
                        f"{self.membrane} and it will be used for solvating the system as a GROMACS triclinic box: \n\t{cryst_info}")
        else:
            self.vectors, self.angles = None, None

        self.cwd = os.getcwd()

    def small_mol_process(self, mol_definition: dict, name: str = "MOL", safe_naming_prefix: str = None):
        """Get parameters for small molecules: ligands, cofactors, ...

        Parameters
        ----------
        mol_definition : dict
            This is a dictionary with:
                * conf -> The path of the small molecule MOL/SDF file [mandatory]. In case that top is provided,
                this must be a .gro, a ValueError will be raised if it is not the case
                the molecule will not get its parameters.

                * top -> GROMACS topology [optional]. Must be a single file topology with all the force field
                information and without the position restraint included, by default None

                * ff:

                    * type -> openff, gaff or espaloma

                    * code -> force field code [optional], by default depending on type

                        * openff -> openff_unconstrained-2.0.0.offxml

                        * gaff -> gaff-2.11

                        * espaloma -> espaloma-0.3.1


        name : str, optional
            Name to give, by default "MOL"
        safe_naming_prefix : str, optional
            This is used to be sure that there will not happen any naming conflict in hte topologies, by default None
        name and safe_naming_prefix will only be used if top is not provided in mol_definition.

        Returns
        -------
        object
            The BioSimSpace system
        """
        force_field_code_default = {
            'openff': 'openff_unconstrained-2.0.0.offxml',
            'gaff': 'gaff-2.11',
            'espaloma': 'espaloma-0.3.1'
        }
        dict_to_work = {
            'top': None,
            'ff': {
                'type': 'openff',
                'code': None,
            }
        }
        if mol_definition:
            recursive_update_dict(dict_to_work, mol_definition)
            dict_to_work['ff']['type'] = str(dict_to_work['ff']['type']).lower()
            if dict_to_work['ff']['type'] not in force_field_code_default and not dict_to_work['top']:
                raise ValueError(f"Molecule {dict_to_work} has non valid type for the force field. Choose from {force_field_code_default.keys()}.")
            # Plug back the default option in case that the user defined None for the code but type was provided correctly
            if not dict_to_work['ff']['code']:
                dict_to_work['ff']['code'] = force_field_code_default[dict_to_work['ff']['type']]
        else:
            raise ValueError(f"Molecule {mol_definition} has a wrong definition.")
        if dict_to_work['conf']:
            if dict_to_work['top']:
                logger.info(f"Using supplied: {dict_to_work['top']} for {dict_to_work['conf']}")
            else:
                logger.info(f"Getting {dict_to_work['ff']['code']} (type = {dict_to_work['ff']['type']}) parameters for: {dict_to_work['conf']}")
        else:
            raise ValueError(f"Molecule {mol_definition} has a wrong configuration")
        # Set flag to False by default
        provided_top_flag = False
        if dict_to_work['top']:
            top_file = os.path.abspath(dict_to_work['top'])
            # In case the user provided a top, set the flag to True
            provided_top_flag = True
            if os.path.splitext(dict_to_work['conf'])[-1] == '.gro':
                gro_file = os.path.abspath(dict_to_work['conf'])
            else:
                raise ValueError("For safety reasons, if top is provided for small molecule; "
                                 f"the gro file must be provided. Provided: {dict_to_work['conf']}.")
        else:
            parameterizer = Parameterize(
                force_field_code=dict_to_work['ff']['code'],
                force_field_type=dict_to_work['ff']['type'],
                ext_types=['top', 'gro'],
                hmr_factor=self.hmr_factor,
                overwrite=True,
                safe_naming_prefix=safe_naming_prefix,
                out_dir=self.wd,
            )
            # Actually you can pass to parameterize Chem.rdchem.Mol, *.inchi, *.smi, *.mol, *.mol2, *.sdf
            parameterizer(input_mol=dict_to_work['conf'], mol_resi_name=name)

            top_file = os.path.join(self.wd, f"{name}.top")
            gro_file = os.path.join(self.wd, f"{name}.gro")

        parmed_system = readParmEDMolecule(top_file=top_file, gro_file=gro_file)
        if provided_top_flag and self.hmr_factor:
            HMassRepartition(parmed_system, self.hmr_factor).execute()
        return parmed_system

    def gmx_process(self, mol_definition: dict, is_membrane: bool = False):
        """Used to process the biomolecules compatibles.
        By default it will use amber99sb-ildn (protein, DNA, ..) Slipids_2020 (membrane).
        However, these setups are overwrite by the definitions on mol_definition 
        Parameters
        ----------
        mol_definition : dict
            This dictionary is wither self.protein or self.membrane, its definition is accessible at the 
            constructor-method's documentation.
        is_membrane : bool, optional
            If True, Slipids_2020 will be set as internal default instead of amber99sb-ildn, by default False
        Returns
        -------
        object
            A par
        """
        # Setting default parameters
        dict_to_work = {
            'top': None,
            'ff': {
                'code': 'Slipids_2020' if is_membrane else 'amber99sb-ildn',
            }
        }
        if mol_definition:
            recursive_update_dict(dict_to_work, mol_definition)
        else:
            return None
        if dict_to_work['conf']:
            dict_to_work['conf'] = os.path.abspath(dict_to_work['conf'])
            name, ext = os.path.splitext(os.path.basename(dict_to_work['conf']))
            if dict_to_work['top']:
                # Convert to absolute paths
                dict_to_work['top'] = os.path.abspath(dict_to_work['top'])
                logger.info(f"Using supplied: {dict_to_work['top']} for {dict_to_work['conf']}")
            else:
                logger.info(f"Getting {dict_to_work['ff']['code']} parameters for: {dict_to_work['conf']}")

        else:
            return None

        gro_out = os.path.join(self.wd, f'{name}.gro')
        top_out = os.path.join(self.wd, f'{name}.top')
        posre_out = os.path.join(self.wd, f'{name}_posre.itp')

        if is_membrane:
            if dict_to_work['ff']['code'] == 'Slipids_2020':
                # Retrieve internal Slipids_2020 only if the user did not provided this force field
                if self.custom_ff_path:
                    if 'Slipids_2020' not in os.listdir(self.custom_ff_path):
                        get_gmx_ff('Slipids_2020', out_dir=self.wd)
                else:
                    get_gmx_ff('Slipids_2020', out_dir=self.wd)

        # TODO Something strange is going on with the posre files. Those are not been used, however, the include statement should be in the topology
        # and becasue I call fix_topology on __call__ the include section of the position
        # restraint should be duplicated but it is not the case, not sure why
        os.chdir(self.wd)
        if dict_to_work['top']:
            shutil.copy(dict_to_work['top'], top_out)
            if ext == '.pdb':
                run(f"gmx editconf -f {dict_to_work['conf']} -o {gro_out}")
            elif ext == '.gro':
                shutil.copy(dict_to_work['conf'], gro_out)
            else:
                raise ValueError(f"Extension of {dict_to_work['conf']} must be .gro or .pdb")
        else:
            if is_membrane:
                run(f"gmx pdb2gmx -f {dict_to_work['conf']} -ff {dict_to_work['ff']['code']} -water none -o {gro_out} -p {top_out} -i {posre_out}")
            else:
                env_prefix = os.environ["CONDA_PREFIX"]
                fixed_pdb = os.path.join(self.wd, f"{name}_fixed.pdb")
                run(f"{env_prefix}/bin/pdbfixer {dict_to_work['conf']} --output={fixed_pdb} --add-atoms=all --replace-nonstandard")
                run(f"gmx pdb2gmx -f {fixed_pdb} -merge all -ff {dict_to_work['ff']['code']} "
                    f"-water none -o {gro_out} -p {top_out} -i {posre_out} -ignh")
        os.chdir(self.cwd)

        # TODO and readParmEDMolecule fails with amber99sb-start-ildn
        system = readParmEDMolecule(top_file=top_out, gro_file=gro_out)
        if self.hmr_factor:
            HMassRepartition(system, self.hmr_factor).execute()
        system.write(os.path.join(self.wd, f'{name}_final.top'))

        return system

    def make_system(self, ligand_definition: dict):
        """Create self.sys_ligand, self.sys_cofactor, self.sys_protein, self.sys_membrane
        and self.md_system (the combination of the available components). In case
        that the class was already called, it will be assumed that self.sys_cofactor, self.sys_protein, self.sys_membrane
        ere already calculated, only self.sys_ligand will be updated as well self.md_system

        Parameters
        ----------
        ligand_definition : dict
            This is a dictionary with. Its definition is the same as mol_definition of the methods self.small_mol_process.
        """
        logger.info("Processing system components")
        self.sys_ligand = self.small_mol_process(
            mol_definition=ligand_definition,
            name="LIG",
            safe_naming_prefix='x')

        # Only if the class has not yet called the full build will be carry out.
        if self.__self_was_called:
            logger.info("Reusing components from cache")
        else:
            if self.cofactor:
                self.sys_cofactor = self.small_mol_process(
                    mol_definition=self.cofactor,
                    name="COF",
                    safe_naming_prefix='z')
            else:
                self.sys_cofactor = None
            self.sys_protein = self.gmx_process(mol_definition=self.protein)
            self.sys_membrane = self.gmx_process(mol_definition=self.membrane, is_membrane=True)
        logger.info("Merging Components")
        # Cofactor at the end in case is a water molecule, not complains from GROMACS after solvation
        self.md_system = system_combiner(protein=self.sys_protein, membrane=self.sys_membrane, ligand=self.sys_ligand, cofactor=self.sys_cofactor)

    def clean(self):
        """Small cleaner, the intermediate steps saved on builder_dir will be deleted
        """
        os.chdir(self.cwd)

        try:
            shutil.rmtree(self.wd)
        except FileNotFoundError:
            pass

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, exception_traceback):
        self.clean()

    def __call__(self, ligand_definition: Union[dict, PathLike], out_dir: str = 'abfe'):
        """The call implementation. It identify if it is needed to build
        all the components of the systems,
        In case that the class was already called, it will assume that all the components of the system,
        with the exception of the ligand, were already builded. This is useful to call the class
        on several ligands that share the same components: protein, membrane and cofactor

        Parameters
        ----------
        ligand_definition : Union[dict, PathLike]
            In case of dictionary, it should have:

                * conf -> The path of the small molecule MOL/SDF file [mandatory]. In case that top is provided,
                this must be a .gro, a ValueError will be raised if it is not the case
                the molecule will not get its parameters.

                * top -> GROMACS topology [optional]. Must be a single file topology with all the force field
                information and without the position restraint included, by default None

                * ff:

                    * type -> openff, gaff or espaloma

                    * code -> force field code [optional], by default depending on type

                        * openff -> openff_unconstrained-2.0.0.offxml

                        * gaff -> gaff-2.11

                        * espaloma -> espaloma-0.3.1

                In case of PathLike:

                * The path of the small MOL/SDF molecule file

        out_dir : str, optional
            Where you would like to export the generated files, by default 'abfe'
        """
        logger.info(39*"-")
        if not isinstance(ligand_definition, dict):
            ligand_definition = {
                'conf': ligand_definition
            }
        logger.info(f"Processing ligand: {ligand_definition['conf']}")
        # Update (on multiple calls) or just create the out_dir (first call)
        self.out_dir = out_dir
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)

        # Construct MD system:
        self.make_system(ligand_definition)
        system_dir = os.path.join(self.wd, 'system')
        ligand_dir = os.path.join(self.wd, 'ligand')

        logger.info(f"Solvating with {self.water_model}:")
        if self.membrane:
            # TODO: this is the easiest way to implement the position restraints changing the restraints
            # during different steps. However, we are not using different restraints for different molecules
            # what might be needed for some systems. That will take some coding in order to identify the molecules
            f_xyz_complex = 3*['POSRES_DYNAMIC']
        else:
            f_xyz_complex = 3*[2500]

        with solvent.Solvate(self.water_model, builder_dir=os.path.join(self.wd, '.solvating')) as SolObj:

            logger.info(f"Ligand in: {ligand_dir}")
            SolObj(structure=self.sys_ligand, bt='octahedron', d=1.5, out_dir=ligand_dir, out_name='solvated', f_xyz=3*[2500])

            logger.info(f"Complex in: {system_dir}")
            settles_to_constraints_on = None
            if self.cofactor:
                if 'is_water' in self.cofactor:
                    if self.cofactor['is_water']:
                        warnings.warn(f'Provided cofactor {self.cofactor} was labeled as water (is_water = True). '
                                      'So, its settles section (if any), will be changed to tip3p-like triangular constraints. '
                                      'Check here for more information: '
                                      'https://gromacs.bioexcel.eu/t/how-to-treat-specific-water-molecules-as-ligand/3470/9')
                        settles_to_constraints_on = 'COF'

            if self.membrane:
                SolObj(structure=self.md_system, bt='triclinic', box=self.vectors, angles=self.angles,
                       out_dir=system_dir, out_name='solvated', f_xyz=f_xyz_complex, settles_to_constraints_on=settles_to_constraints_on)
            else:
                SolObj(structure=self.md_system, bt='octahedron', d=1.5, out_dir=system_dir, out_name='solvated',
                       f_xyz=f_xyz_complex, settles_to_constraints_on=settles_to_constraints_on)

        # Make index file in case of membrane systems
        if self.membrane:
            solvent.index_for_membrane_system(
                configuration_file=os.path.join(system_dir, "solvated.gro"),
                ndxout=os.path.join(system_dir, "index.ndx"),
                lignad_name='LIG',
                cofactor_name='COF' if self.cofactor else None,
                cofactor_on_protein=self.cofactor_on_protein,
            )
        else:
            # Create a dummy index.ndx file. It is needed for the Snakemake workflow
            open(os.path.join(system_dir, "index.ndx"), "w").close()

        # Construct ABFE system:
        logger.info(f"Final build of ABFE directory on: {self.out_dir}")
        make_abfe_dir(out_dir=self.out_dir, ligand_dir=ligand_dir, sys_dir=system_dir)

        # Change state
        self.__self_was_called = True
        logger.info("--------- Building Completed ----------\n")

#############################################################################################


if __name__ == "__main__":
    pass
