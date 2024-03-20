#!/usr/bin/env python
import os
import subprocess
import tempfile
from pathlib import Path
from typing import List, Union, Tuple

PathLike = Union[os.PathLike, str, bytes]


class DotDict:
    """A simple implementation of dot-access dict"""
    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            if isinstance(value, dict):
                self.__dict__[key] = DotDict(**value)
            else:
                self.__dict__[key] = value

    def __repr__(self) -> str:
        return str(self.__dict__)

    def to_dict(self):
        return self.__dict__


def run(command: str, shell: bool = True, executable: str = '/bin/bash', interactive: bool = False) -> subprocess.CompletedProcess:
    """A simple wrapper around subprocess.Popen/subprocess.run

    Parameters
    ----------
    command : str
        The command line to be executed
    shell : bool, optional
        Create a shell section, by default True
    executable : str, optional
        what executable to use, pass `sys.executable` to check yours, by default '/bin/bash'
    interactive : bool, optional
        To interact with the command, by default False. If True, you can access stdout and stderr of the returned process.

    Returns
    -------
    subprocess.CompletedProcess
        The process

    Raises
    ------
    RuntimeError
        In case that the command fails, the error is raised in a nice way
    """
    if interactive:
        process = subprocess.run(command, shell=shell, executable=executable)
        returncode = process.returncode
        if returncode != 0:
            raise RuntimeError(f'Command {command} returned non-zero exit status {returncode}')
    else:
        process = subprocess.run(command, shell=shell, executable=executable, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        returncode = process.returncode

        if returncode != 0:
            print(f'Command {command} returned non-zero exit status {returncode}')
            raise RuntimeError(process.stderr)

    return process


def gmx_command(load_dependencies: List[str] = None, interactive: bool = False, stdout_file: PathLike = None):
    """Lazy wrapper of gmx commands

    Parameters
    ----------
    load_dependencies : List[str]
        It is used in case some previous loading steps are needed;
        e.g: ['source /groups/CBG/opt/spack-0.18.1/shared.bash', 'module load sandybridge/gromacs/2022.4']
    interactive : bool
        In case, and interactive section is desired, by default False
    stdout_file : bool
        IF provided, it will append to the command ` >& {stdout_file}`, by default None


    A typical function will be:

    Example
    -------
    .. ipython:: python

        from abfe.utils import tools
        @tools.gmx_command()
        def mdrun(
                **kwargs
                ):...

    The important parts are:

    #. The name of the function must be the name of the gmx command, for example mdrun, grompp, etc.
    #. You must return the local variables of the function
    #. The names of the keywords are exactly the same name as got it by the respective function.
    #. For flags, a boolean will be provided as value, for example v = True, if you want to be verbose.
    """
    def decorator(gmx_function: object):
        def wrapper(**kwargs):
            if load_dependencies:
                cmd = " && ".join(load_dependencies)
                cmd += " && "
            else:
                cmd = ''
            cmd += f"gmx {gmx_function.__name__}"
            for key in kwargs:
                value = kwargs[key]
                if value:
                    if isinstance(value, bool):
                        cmd += f" -{key}"
                    else:
                        cmd += f" -{key} {value}"
            if stdout_file:
                cmd += f" >& {stdout_file}"
                if interactive:
                    raise RuntimeError("stdout_file argument is not compatible with interactive flag")
            if interactive:
                return run(cmd, interactive=True)
            else:
                return run(cmd)
        return wrapper
    return decorator


def gmx_runner(mdp: PathLike, topology: PathLike, structure: PathLike, checkpoint: PathLike = None, index: PathLike = None,
               nthreads: int = 12, load_dependencies: List[str] = None, run_dir: PathLike = '.', **mdrun_extra):
    """This function create the tpr file based on the input provided
    And run the simulation.
    Note: During the tpr creation maxwarn = 2 (TODO: remove it in the future)

    The following commands will be executed by default:

    gmx grompp -f {mdp} -c {structure} -r {structure} -p {topology} -o {mdp-name}.tpr -maxwarn 2
    gmx mdrun -nt 12 -deffnm {mdp-name}

    mdrun will update the command based on mdrun_extra. You can also suppress the use of nt and/or deffnm passing them as
    False and construct your own mdrun command:
    e.g. gmx_runner(mdp='emin.mdp', topology='ligand.top',structure='ligand.gro', deffnm = False, cpi = True, s = 'emin.tpr', o = 'emin2', c = 'emin3')

    The last will give:
    gmx mdrun -nt 12 -cpi -s emin.tpr -o emi666 -c emi55

    Parameters
    ----------
    mdp : str
        The path to the MDP file. The name of the file will be used for the tpr and for the files generated during mdrun.
    topology : PathLike
        GMX topology file
    structure : PathLike
        The PDB, GRO, etc structure of the system
    checkpoint : PathLike
        Full precision trajectory: trr cpt tng, by default None> if given will be used on grompp with the flag `-t {checkpoint}`
    index : PathLike
        A GMX index to be used on grompp, by default None
    nthreads : int, optional
        Number of threads to run, by default 12
    load_dependencies : List[str], optional
        It is used in case some previous loading steps are needed;
        e.g: ['source /groups/CBG/opt/spack-0.18.1/shared.bash', 'module load sandybridge/gromacs/2022.4'], by default None
    run_dir : PathLike, optional
        Where the simulation should run (write files). If it does not exist will be created, by default '.'
    **mdrun_extra : any
        Any valid keyword for mdrun. flags are passing as boolean. E.g: cpi = True
    """
    # Create run directory on demand
    makedirs(run_dir)

    name = os.path.splitext(os.path.basename(mdp))[0]

    # Because of how snakemake handles environmental variables that
    # are used by GROMACS (https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html)
    # We have to hard code the unset of some of them
    # TODO, check the implications of such modifications. I hope that only affects the specific rule where GROMACS is called
    hard_code_dependencies = [
        'unset OMP_NUM_THREADS',
        'unset GOTO_NUM_THREADS',
        'unset OPENBLAS_NUM_THREADS',
        'unset MKL_NUM_THREADS',
        'unset VECLIB_MAXIMUM_THREADS',
        'unset NUMEXPR_NUM_THREADS',
    ]

    @gmx_command(load_dependencies=hard_code_dependencies + load_dependencies)
    def grompp(**kwargs): ...

    @gmx_command(load_dependencies=hard_code_dependencies + load_dependencies, stdout_file=f"{name}.lis")
    def mdrun(**kwargs): ...

    cwd = os.getcwd()
    os.chdir(run_dir)

    grompp_extra = {}
    if checkpoint:
        grompp_extra['t'] = checkpoint
    if index:
        grompp_extra['n'] = index

    # TODO, I do not like to use the maxwarn keyword hardcoded.
    grompp(f=f"{mdp}", c=structure, r=structure, p=topology, o=f"{name}.tpr", maxwarn=2, **grompp_extra)

    mdrun_kwargs = {
        # TODO DEBUG
        # "ntomp": nthreads,
        "nt": nthreads,
        # TODO: this flag is deprecate in new GROMACS versions and in 2024 is not longer available
        # This means that I have to build the mdrun command at rule level.
        "deffnm": name,
    }
    if mdrun_extra:
        mdrun_kwargs.update(mdrun_extra)
    mdrun(**mdrun_kwargs)

    os.chdir(cwd)


def paths_exist(paths: List, raise_error: bool = False, out: Union[str, None] = None) -> None:
    """Check that the paths exist

    Parameters
    ----------
    paths : List
        A list of paths
    raise_error : bool, optional
        If True will raise a RuntimeError when any path doe not exist, by default False
    out : Union[str, None], optional
        In case that all files exist and out is st to some file; the existence of this file could be
        used as a check that all paths exist (useful for sanekemake), by default None

    Raises
    ------
    RuntimeError
        In case some path does not exist and rasie_error = True
    """
    check = True
    for path in paths:
        if not os.path.exists(path):
            check = False
            msg = f"Missing path/file: {path}"
            if raise_error:
                raise RuntimeError(msg)
            else:
                print(msg)
    if out and check:
        open(out, "w").close()


def list_if_dir(path: PathLike = '.'):
    return [item for item in os.listdir(path) if os.path.isdir(os.path.join(path, item))]


def list_if_file(path: PathLike = '.', ext: str = None) -> List[str]:
    """Dir all the files in path

    Parameters
    ----------
    path : PathLike, optional
        Path to look for the file, by default '.'
    ext : str, optional
        Th extension of the file, for example: "py", "sh", "txt", by default None

    Returns
    -------
    List[str]
        The list of file names
    """
    files = [item for item in os.listdir(path) if os.path.isfile(os.path.join(path, item))]
    if ext:
        files = [file for file in files if os.path.splitext(file)[-1] == f".{ext}"]
    return files


def is_file_inside_directory(directory_path, file_path):
    # Convert paths to Path objects
    directory_path = Path(directory_path).resolve()
    file_path = Path(file_path).resolve()

    # Check if the file path starts with the directory path
    # print(file_path.parts)
    return file_path.parts[:len(directory_path.parts)] == directory_path.parts


def find_xtc(root_path: PathLike, exclude_suffixes: List[str] = None) -> List[PathLike]:
    """Find all the files with the extension .xtc that does not have any
    parent directory with any exclude_suffixes. If the name if the xtc file
    has as suffix some of the ones specified in excluded_suffixes, it will also
    discarded as well.

    Parameters
    ----------
    root_path : PathLike
        Root path to look for XTC files
    exclude_suffixes : List[str], optional
        list of suffixes to exclude from wither parent directories or the XTC files themself
        , by default None

    Returns
    -------
    List[PathLike]
        LIst of XTC file paths
    """
    xtc_files = Path(root_path).resolve().rglob('*.xtc')
    if exclude_suffixes:
        exclude_suffixes = tuple(exclude_suffixes)
        xtc_files_filtered = []
        for xtc_file in xtc_files:
            components = [xtc_file] + list(xtc_file.parents)
            # any parent directories or the file itself has any of exclude_suffixes
            test = any([True if component.endswith(exclude_suffixes) else False for component in components])
            if not test:
                xtc_files_filtered.append(xtc_file)
        return xtc_files_filtered
    else:
        return xtc_files


def archive(root_path: PathLike, exclude_suffixes: List[str] = None, name: str = 'archive',
            compress_type: str = 'gz', remove_dirs: bool = False, out_check_file: bool = True):
    """Recursively archive root_path. Directories and/or files with ny suffixes from
     exclude_suffixes are ignored . It creates a tar file with the XTC files (without compress)
    and a main_project.tar.{compress_type} with the rest of directories. Compression will only be applied to
    those files included in main_project.tar.{compress_type}. In-house benchmark showed a compress
    rate close to for a abfe campaign 1.8 using gz compression
    (data taken from MCL1).
    139 GB to 77 GB

    .. warning::
        It may be that the function fail because the directory is too large, in this case you must split the directory,
        this was the case for the p38 campaign (https://github.com/openforcefield/protein-ligand-benchmark) with 3 replicas

        BE AWARE OF THE IMPLICATION TO DELETE A SIMULATION DIRECTORY with the option remove_dirs = True

    In-house benchmark showed:

    +-------+-------+-------+
    |       | Time  | Space |
    +=======+=======+=======+
    | bz2   | x s   | x MB  |
    +-------+-------+-------+
    | gz    | x s   | x MB  |
    +-------+-------+-------+
    | xz    | x s   | x MB  |
    +-------+-------+-------+

    Parameters
    ----------
    root_path : PathLike
        The root path for which all the dirs will be compressed
    exclude_suffixes : List[PathLike], optional
        List of suffix to exclude for compression either directories or files. The endswith method will be applied
        Use case example could be: [.snakemake, .log, .edr, .lis, .err]. In this case the directory .snakemake
        will be ignored and all the files with the specified extensions.
    name : str, optional
        Output name of the archive file, by default 'archive'
    compress_type : str, optional
        Type of compression to use, tar, gz, bz2 and xz are possible, by default 'gz'
    remove_dirs : bool, optional
        Remove compressed root_path, by default False
    out_check_file : bool, optional
        If the archive worked as expected, a file {name}_safe_remove.check
        will be written, by default True

    Raises
    ------
    FileNotFoundError
        If root_path does not exist
    ValueError
        Incorrect compress_type
    ValueError
        If the provided name lays on root_path, this is not expected.
    """
    import shutil
    import tarfile

    # Ensure the provided path exists
    if not os.path.exists(root_path):
        raise FileNotFoundError(f"Directory '{root_path}' does not exist.")

    # Define the name of the compressed file
    compress_type = compress_type.lower()
    valid_tar_exts = ['tar', 'gz', 'bz2', 'xz']
    if compress_type not in valid_tar_exts:
        raise ValueError(f"Unsupported compression type ({compress_type}). Use: {' '.join(valid_tar_exts)}.")
    if is_file_inside_directory(root_path, f"{name}.tar"):
        raise ValueError(f"Invalid {name=}. It lays in {root_path=}")

    # Convert to list
    if exclude_suffixes:
        exclude_suffixes = list(set(exclude_suffixes))
    else:
        exclude_suffixes = []

    # Find and create a separate archive for XTC files
    xtc_files = find_xtc(root_path=root_path, exclude_suffixes=exclude_suffixes)

    with tarfile.open(f"{name}.tar", 'w:tar') as project_archive:
        if xtc_files:
            for xtc_file in xtc_files:
                xtc_file_path = os.path.join(root_path, xtc_file)
                print(f"Adding XTC: {xtc_file}")
                project_archive.add(xtc_file_path, arcname=xtc_file)

        with tempfile.TemporaryDirectory(prefix='.main_archive', dir='.') as tmpdir:
            # Create the compressed main_archive
            arcname = 'main_project.tar'
            if compress_type != 'tar':
                arcname += f".{compress_type}"
            # For the main archive always exclude XTC files
            internal_excluded_ext = tuple(set(exclude_suffixes + ['.xtc']))
            with tarfile.open(os.path.join(tmpdir, arcname), f'w:{compress_type}') as main_archive:
                for root, _, files in os.walk(root_path):
                    for file in files:
                        if not file.endswith(internal_excluded_ext):
                            file_path = os.path.relpath(os.path.join(root, file), root_path)
                            main_archive.add(os.path.join(root, file), arcname=file_path)
            print(f"Adding: {arcname}")
            project_archive.add(os.path.join(tmpdir, arcname), arcname=arcname)

    if out_check_file:
        with open(f"{name}_safe_remove.check", "w") as f:
            f.write('All files were successfully archived!')
    # Optionally remove the source directories
    if remove_dirs:
        # TODO Check that f"{name}.tar" does not lay in root_dir
        print("Cleaning after compression:")
        print(f"Removing: {root_path}")
        shutil.rmtree(root_path)


def _filter_helper(TarInfo: str, suffix: Tuple[str], prefix: Tuple[str] = ('main_project.tar')):
    if suffix:
        if TarInfo.name.endswith(suffix):
            return TarInfo
        else:
            if prefix:
                if TarInfo.name.startswith(prefix):
                    return TarInfo
                else:
                    return None
            return None
    return TarInfo


def unarchive(archive_file: PathLike, target_path: PathLike,
              only_with_suffix: Union[None, List[str]] = None, prefix: Tuple[str] = ('main_project.tar')):
    """It unarchive a project archived by the function :meth:`abfe.utils.tools.archive`

    Parameters
    ----------
    archive_file : PathLike
        Archived project
    target_path : PathLike
        Out path to unarchive
    only_with_suffix : Union[None, List[str]]
        Only extract those files that present the suffix
    """
    import tarfile

    # Ensure the target directory exists
    target_path = os.path.abspath(target_path)
    os.makedirs(target_path, exist_ok=True)

    # Convert to list
    if only_with_suffix:
        # Addint the main_project.tar
        only_with_suffix = tuple(only_with_suffix)
    else:
        only_with_suffix = tuple()

    # Create a temporary directory for extracting the main compressed archive
    with tempfile.TemporaryDirectory(prefix='.unarchive_main', dir='.') as tmpdir:
        # Extract the XTC archive first
        with tarfile.open(archive_file, 'r') as archive:
            for member in archive.getmembers():
                member = _filter_helper(member, only_with_suffix, prefix=prefix)
                if member:
                    print(f"Decompressing: {member.name}")
                    if member.name.startswith('main_project.tar'):
                        compress_type = member.name.split('.')[-1]
                        main_archive_path = os.path.join(tmpdir, member.name)
                        # TODO
                        # This step is extremely painful, at some point you have
                        # twice the size of the original archive file
                        # The solution is to have separately for xtc and main project, but then we have two files
                        # not soe "archive", but either we solve the clean archive, or improve the unarchive.
                        archive.extract(member, tmpdir)
                        with tarfile.open(main_archive_path, f'r:{compress_type}') as main_archive:
                            for main_member in main_archive.getmembers():
                                main_member = _filter_helper(main_member, only_with_suffix, prefix=prefix)
                                if main_member:
                                    print(f"Decompressing: {main_member.name}")
                                    main_archive.extract(main_member, target_path)
                    else:
                        archive.extract(member, target_path)


def makedirs(path):
    if not os.path.exists(path):
        os.makedirs(path)


def recursive_update_dict(original_dict: dict, update_dict: dict) -> None:
    for key, value in update_dict.items():
        if isinstance(value, dict) and key in original_dict and isinstance(original_dict[key], dict):
            recursive_update_dict(original_dict[key], value)
        else:
            original_dict[key] = value


def config_validator(global_config: dict) -> List:
    """It checks for the validity of the global config.
    This dictionary is usually passed to :meth:`abfe.calculate_abfe.calculate_abfe`

    Parameters
    ----------
    global_config : dict
        The configuration of the ABFE workflow

    Returns
    -------
    List[bool,str]
        result[0], True if pass all the checks. False otherwise.
        result[1], Extra information.
    """

    # Checking cluster
    if 'cluster' not in global_config:
        return False, "Cluster configuration is missing"

    if 'type' not in global_config['cluster']:
        return False, "Cluster type is missing"

    if 'options' not in global_config['cluster']:
        return False, "Cluster configuration is valid, but no cluster options provided"

    if 'calculation' not in global_config['cluster']['options']:
        return False, "cluster/options configuration is valid, but no cluster/options/calculation provided"

    # Setting up default extra mdrun and job dependencies in case it was not provided
    if "extra_directives" in global_config:
        if "dependencies" not in global_config["extra_directives"]:
            global_config["extra_directives"]["dependencies"] = []
        if "mdrun" not in global_config["extra_directives"]:
            global_config["extra_directives"]["mdrun"] = {
                'ligand': {},
                'complex': {},
                'all': {}
            }
    else:
        global_config["extra_directives"] = {
            "dependencies": [],
            "mdrun": {
                'ligand': {},
                'complex': {},
                'all': {}
            },
        }
    valid_mdrun = ["ligand", "complex", "all"]
    # In case that 'extra_directives/mdrun/key' was not defined
    for key in valid_mdrun:
        if key not in global_config["extra_directives"]["mdrun"]:
            global_config["extra_directives"]["mdrun"][key] = {}

    # Check that mdrun is valid
    valid_mdrun = ["ligand", "complex", "all"]
    for key in global_config["extra_directives"]["mdrun"]:
        if key not in valid_mdrun:
            return False, f"extra_directives/mdrun/{key} is not valid, you must select one of valid mdrun options {valid_mdrun}"

        # Here we use as base keywords the one defined in all
        # And then, for ligand and complex, update those based on th user input
        # In other words, ligand and complex will use the all definition updated by their own keywords.
        if key != 'all':
            key_all = global_config["extra_directives"]["mdrun"]['all'].copy()
            key_all.update(global_config["extra_directives"]["mdrun"][key])
            global_config["extra_directives"]["mdrun"][key] = key_all

        # Always allow continuation in case the user did not defined
        if "cpi" not in global_config["extra_directives"]["mdrun"][key]:
            global_config["extra_directives"]["mdrun"][key]['cpi'] = True

    # After the update keywords, keep all is not needed any more
    del global_config["extra_directives"]["mdrun"]['all']

    return True, "Cluster configuration is valid"


if __name__ == "__main__":
    pass
