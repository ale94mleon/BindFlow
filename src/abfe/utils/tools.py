#!/usr/bin/env python
import subprocess, os
from typing import List, Union

PathLike = Union[os.PathLike, str, bytes]

def run(command:str, shell:bool = True, executable:str = '/bin/bash', interactive:bool = False) -> subprocess.CompletedProcess:
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
        process = subprocess.run(command, shell = shell, executable = executable, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text = True)
        returncode = process.returncode
    
        if returncode != 0:
            print(f'Command {command} returned non-zero exit status {returncode}')
            raise RuntimeError(process.stderr)

    return process

def gmx_command(load_dependencies:List[str] = None, interactive:bool = False):
    """Lazy wrapper of gmx commands

    Parameters
    ----------
    load_dependencies : List[str]
        It is used in case some previous loading steps are needed;
        e.g: ['source /groups/CBG/opt/spack-0.18.1/shared.bash', 'module load sandybridge/gromacs/2022.4']
    A typical function will be:
    
    Example
    -------
    .. ipython:: python

        def mdrun(
                **kwargs
                ):...
    
    The important parts are:

    #. The name of the function must be the name of the gmx command, for example mdrun, grompp, etc.
    #. You must return the local variables of the function
    #. The names of the keywords are exactly the same name as got it by the respective function.
    #. For flags, a boolean will be provided as value, for example v = True, if you want to be verbose.
    """   
    def decorator(gmx_function:object):
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
            print(cmd)
            if interactive:
                return run(cmd, interactive=True)
            else:
                return run(cmd)
        return wrapper
    return decorator

def gmx_runner(mdp:PathLike, topology:PathLike, structure:PathLike, checkpoint:PathLike = None, nthreads:int = 12, load_dependencies:List[str] = None, run_dir:PathLike = '.', **mdrun_extra):
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
    nthreads : int, optional
        Number of threads to run, by default 12
    load_dependencies : List[str], optional
        It is used in case some previous loading steps are needed;
        e.g: ['source /groups/CBG/opt/spack-0.18.1/shared.bash', 'module load sandybridge/gromacs/2022.4'], by default None
    run_dir : PathLike, optional
        Where the simulation should run (write files). If it doe snot exist will be created, by default '.'
    **mdrun_extra : any
        Any valid keyword for mdrun. flags are passing as boolean. E.g: cpi = True
    """
    # Create run directory on demand
    makedirs(run_dir)

    @gmx_command(load_dependencies=load_dependencies)
    def grompp(**kwargs):...
    
    @gmx_command(load_dependencies=load_dependencies)
    def mdrun(**kwargs):...
    
    name = os.path.splitext(os.path.basename(mdp))[0]
    cwd = os.getcwd()
    os.chdir(run_dir)
    # TODO, I do not like to use the maxwarn keyword hardcoded.
    
    if checkpoint:
        grompp_extra = {'t':checkpoint}
    else:
        grompp_extra = {}
    grompp(f = f"{mdp}", c = structure, r = structure, p = topology, o = f"{name}.tpr", maxwarn = 2, **grompp_extra)
    
    mdrun_kwargs = {
        "nt": nthreads,
        "deffnm": name,
    }
    if mdrun_extra:
        mdrun_kwargs.update(mdrun_extra)
    mdrun(**mdrun_kwargs)

    os.chdir(cwd)

def list_if_dir(path = '.'):
    return [item for item in os.listdir(path) if os.path.isdir(os.path.join(path, item))]

def list_if_file(path:PathLike = '.', ext:str = None) -> List[str]:
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

def makedirs(path):
    if os.path.exists(path):
        pass
    else:
        os.makedirs(path,exist_ok=True)

def config_validator(global_config:dict) -> List:
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
        if not "dependencies" in global_config["extra_directives"]:
            global_config["extra_directives"]["dependencies"] = []
        if not "mdrun" in global_config["extra_directives"]:
            global_config["extra_directives"]["mdrun"] = {}    
    else:
        global_config["extra_directives"] = {
            "dependencies": [],
            "mdrun": {},
        }


    return True, "Cluster configuration is valid"



if __name__ == "__main__":
    pass
    glolba_config = {
        "cluster": {
            # "type": "slurm",
            "options":{
                "calculation": {
                    1:1
                }
            }
        }
    }
    print(config_validator(global_config=glolba_config))