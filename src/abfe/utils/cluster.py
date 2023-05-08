#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from typing import List
import json, os, subprocess

_SBATCH_KEYWORDS = {
    #https://slurm.schedmd.com/sbatch.html
    'A':'account','account':'account',
    'acctg-freq':'acctg-freq',
    'a':'array','array':'array',
    'batch':'batch',
    'bb':'bb',
    'bbf':'bbf',
    'b':'begin', 'begin': 'begin',
    'D':'chdir', 'chdir':'chdir',
    'cluster-constraint':'cluster-constraint',
    'M':'clusters', 'clusters':'clusters',
    'comment':'comment',
    'C':'constraint', 'constraint':'constraint',
    'container':'container',
    'contiguous':'contiguous',
    'S':'core-spec', 'core-spec':'core-spec',
    'cores-per-socket':'cores-per-socket',
    'cpu-freq':'cpu-freq',
    'cpus-per-gpu':'cpus-per-gpu',
    'c':'cpus-per-task', 'cpus-per-task':'cpus-per-task',
    'deadline':'deadline',
    'delay-boot':'delay-boot',
    'd':'dependency', 'dependency':'dependency',
    'm':'distribution', 'distribution':'distribution',
    'e':'error', 'error':'error',
    'x':'exclude', 'exclude':'exclude',
    'exclusive':'exclusive',
    'export':'export',
    'export-file':'export-file',
    'B':'extra-node-info', 'extra-node-info':'extra-node-info',
    'get-user-env':'get-user-env',
    'gid':'gid',
    'gpu-bind':'gpu-bind',
    'gpu-freq':'gpu-freq',
    'G':'gpus', 'gpus':'gpus',
    'gpus-per-node':'gpus-per-node',
    'gpus-per-socket':'gpus-per-socket',
    'gpus-per-task':'gpus-per-task',
    'gres':'gres',
    'gres-flags':'gres-flags',
    'h':'help', 'help':'help',
    'hint':'hint',
    'H':'hold', 'hold':'hold',
    'ignore-pbs':'ignore-pbs',
    'i':'input', 'input':'input',
    'J':'job-name', 'job-name':'job-name',
    'kill-on-invalid-dep':'kill-on-invalid-dep',
    'L':'licenses', 'licenses':'licenses',
    'mail-type':'mail-type',
    'mail-user':'mail-user',
    'mcs-label':'mcs-label',
    'mem':'mem',
    'mem-bind':'mem-bind',
    'mem-per-cpu':'mem-per-cpu',
    'mem-per-gpu':'mem-per-gpu',
    'mincpus':'mincpus',
    'network':'network',
    'nice':'nice',
    'k':'no-kill', 'no-kill':'no-kill',
    'no-requeue':'no-requeue',
    'F':'nodefile', 'nodefile':'nodefile',
    'w':'nodelist','nodelist':'nodelist',
    'N':'nodes', 'nodes':'nodes',
    'n':'ntasks', 'ntasks':'ntasks',
    'ntasks-per-core':'ntasks-per-core',
    'ntasks-per-gpu':'ntasks-per-gpu',
    'ntasks-per-node':'ntasks-per-node',
    'ntasks-per-socket':'ntasks-per-socket',
    'open-mode':'open-mode',
    'o':'output', 'output':'output',
    'O':'overcommit', 'overcommit':'overcommit',
    's':'oversubscribe', 'oversubscribe':'oversubscribe',
    'parsable':'parsable',
    'p':'partition', 'partition':'partition',
    'power':'power',
    'priority':'priority',
    'profile':'profile',
    'propagate':'propagate',
    'q':'qos', 'qos':'qos',
    'Q':'quiet','quiet':'quiet',
    'reboot':'reboot',
    'requeue':'requeue',
    'reservation':'reservation',
    'signal':'signal',
    'sockets-per-node':'sockets-per-node',
    'spread-job':'spread-job',
    'switches':'switches',
    'test-only':'test-only',
    'thread-spec':'thread-spec',
    'threads-per-core':'threads-per-core',
    't':'time', 'time':'time', #Acceptable time formats include "minutes", "minutes:seconds", "hours:minutes:seconds", "days-hours", "days-hours:minutes" and "days-hours:minutes:seconds".
    'tmp':'tmp',
    'uid':'uid',
    'usage':'usage',
    'use-min-nodes':'use-min-nodes',
    'v':'verbose', 'verbose':'verbose',
    'V':'version', 'version':'version',
    'W':'wait', 'wait':'wait',
    'wait-all-nodes':'wait-all-nodes',
    'wckey':'wckey',
    'wrap':'wrap',
}

def remove_initial_dash(string: str) -> str:
    if string.startswith("--"):
        return string[2:]
    elif string.startswith("-"):
        return string[1:]
    else:
        return string

class BaseCluster:
    # Default class variables
    submit_command = "sbatch"
    cancel_command = "scancel"
    config_name = "slurm"
    shebang = "#!/bin/bash"
    job_keyword = "SBATCH"
    
    def __init__(self, scheduler_directives:dict = None, job_extra_directives:List[str] = None, mdrun_extra_directives:dict = None) -> None:
        """ Initializes a new instance of the BaseCluster class.

        Parameters
        ----------
        scheduler_directives :dict, optional
            Dictionary of scheduler directives to be passed to the job script in the form of `scheduler_key:value`, by default None
        job_extra_directives : List[str], optional
            List of extra directives to be added to the job script. This is used if it is needed to launch some modules like spack, conda, etc
            Source files or directories, change directories, by default None
        mdrun_extra_directives : dict, optional
            Dictionary of extra options to be passed to the GROMACS mdrun command in the form of `mdrun_key:value`, by default None
            e.g: {'update':'gpu', 'cpi':True}
        """
        self.scheduler_directives = scheduler_directives
        self.job_extra_directives = job_extra_directives
        # TODO: Control on the GMX options, valid commands
        self._main_mdrun = 'gmx mdrun'
        self.mdrun_extra_directives = mdrun_extra_directives

        self.jobid = None

    
    def set_main_mdrun(self, main_mdrun_cmd:str):
        """This will set that main part of the mdrun command.
        This will be complemented with `self.mdrun_extra_directives` for the final GMX command construction.
        Therefore you must not repeat keywords.

        Parameters
        ----------
        main_mdrun_cmd : str
            The main part of the command; e.g: `gmx mdrun -deffnm production`

        Raises
        ------
        ValueError
            In case that some keywords are already defined on `self.mdrun_extra_directives`
        """
        for item in main_mdrun_cmd.split():
            if f"-{remove_initial_dash(item)}" in self.mdrun_extra_directives:
                raise ValueError(f"{item} is already defined on `mdrun_extra_directives`")
        self._main_mdrun = main_mdrun_cmd

    def __get_full_data(self):
        data = {
            "submit_command": self.__class__.submit_command,
            "cancel_command": self.__class__.cancel_command,
            "config_name": self.__class__.config_name,
            "shebang": self.__class__.shebang,
            "job_keyword": self.__class__.job_keyword,
        }
        data.update(self.__dict__)
        return data
    
    def to_json(self, out_file:str = "cluster.json"):
        """Method to write all the attributes of the BaseCluster class to a JSON file

        Parameters
        ----------
        out_file : str, optional
            Name of the output JSON file, by default "cluster.config".
        """

        with open(out_file, 'w') as f:
            json.dump(self.__get_full_data(), f, indent=4)

    def get_job_script(self) -> str:
        string = self.shebang+"\n"
        
        for scheduler_key in self.scheduler_directives:
            value = self.scheduler_directives[scheduler_key]
            if isinstance(value, bool):
                value = '\n'
            else:
                value = f" {value}\n"
            string +=  f"{self.job_keyword} {scheduler_key}{value}"

        string += "\n"
        
        for item in self.job_extra_directives:
            string += item+"\n"
        string += "\n"

        # Setting the mdrun
        string += self._main_mdrun
        for mdrun_key in self.mdrun_extra_directives:
            value = self.mdrun_extra_directives[mdrun_key]
            if isinstance(value, bool):
                value = ""
            string += f" -{mdrun_key} {value}"
        
        return string
    
    def __repr__(self):
        return f"{self.__class__.__name__}(\n{json.dumps(self.__get_full_data(), indent=5)}\n)"
    
    def submit(self,wd = '.'):
        cwd = os.getcwd()
        os.chdir(wd)
        # Build the job script or the command line of SLURM, in this case is build the slurm command 
        # update jobid as 
        cmd = f"{self.submit_command} "
        os.chdir(cwd)
    
    def get_job_status(self):
        raise NotImplementedError(f"{self.__class__.__name__} this method is not yet implemented")

class SLURMCluster(BaseCluster):
    # Override class variables
    submit_command = "sbatch"
    cancel_command = "scancel"
    config_name = "slurm"
    shebang = "#!/bin/bash"
    job_keyword = "#SBATCH"

    def __init__(self, scheduler_directives:dict = None, state:List[str] = None, mdrun_extra_directives:List[str] = None) -> None:
        """ Initializes a new instance of the SLURMCluster class.

        Parameters
        ----------
        scheduler_directives : List[str], optional
            Dictionary of SLURM directives to be passed to the job script, by default None
        job_extra_directives : List[str], optional
            List of extra directives to be added to the job script. This is used if it is needed to launch some modules like spack, conda, etc
            Source files or directories, change directories, by default None
        mdrun_extra_directives : List[str], optional
            Dictionary of extra options to be passed to the GROMACS mdrun command in the form of `mdrun_key:value`, by default None
    Example
    -------
    .. ipython:: python

        from abfe.utils import cluster
        scheduler_directives = {"partition": "my_partition", "time": "10:00:00", "c": 16, "mem-per-cpu": "4G"}
        job_extra_directives = ["source /groups/CBG/opt/spack-0.18.1/shared.bash", "module load gromacs/2022.4"]
        mdrun_extra_directives = {"bonded": "gpu", "npme": 2, "cpi": True}
        my_cluster = cluster.SLURMCluster(scheduler_directives=scheduler_directives, job_extra_directives=job_extra_directives, mdrun_extra_directives=mdrun_extra_directives)

        # Set the main mdrun command, in this way the class could be used several times for several runs
        my_cluster.set_main_mdrun("gmx-mpi mdrun -deffnm production")

        # Print the class
        print(my_cluster)

        # Print the cluster script
        print(my_cluster.get_job_script())

        # Write the attributes of the object to a JSON file
        my_cluster.to_json("my_cluster.json")
        """
        # Translate scheduler_directives
        translated_scheduler_directives = {}
        for key in scheduler_directives:
            if key not in _SBATCH_KEYWORDS: raise ValueError(f"{key} is not a valid SLURM string key")
            translated_scheduler_directives[_SBATCH_KEYWORDS[key]] = scheduler_directives[key]

        super().__init__(scheduler_directives=translated_scheduler_directives, job_extra_directives=job_extra_directives, mdrun_extra_directives=mdrun_extra_directives)

    def get_job_status(self):
        output = str(subprocess.check_output("sacct -j %s --format State --noheader | head -1 | awk '{print $1}'" % self.jobid, shell=True).strip())

        running_status = ["PENDING", "CONFIGURING", "COMPLETING", "RUNNING", "SUSPENDED"]
        if "COMPLETED" in output:
            status = "COMPLETED"
        elif any(r in output for r in running_status):
            status = "RUNNING"
        else:
             status = "FAILED"
        return status

if __name__ == '__main__':
    pass
    # scheduler_directives = {"partition": "my_partition", "time": "10:00:00", "c": 16, "mem-per-cpu": "4G"}
    # job_extra_directives = ["source /groups/CBG/opt/spack-0.18.1/shared.bash", "module load gromacs/2022.4"]
    # mdrun_extra_directives = {"bonded": "gpu", "npme": 2, "cpi": True}
    # my_cluster = SLURMCluster(scheduler_directives=scheduler_directives, job_extra_directives=job_extra_directives, mdrun_extra_directives=mdrun_extra_directives)
    # print(isinstance(my_cluster, BaseCluster))
    # # Write the attributes of the object to a JSON file
    # # my_cluster.set_main_mdrun("gmx-mpi mdrun -deffnm production")
    # # print(my_cluster)
    # # print(my_cluster.get_job_script())
