import json
import os
import stat
from abfe.utils.cluster import _SBATCH_KEYWORDS
from abfe.utils import tools
from abfe.utils.tools import PathLike
from abc import ABC, abstractmethod

class Scheduler(ABC):
    # Default class variables
    submit_command = None
    cancel_command = None
    shebang = None
    job_keyword = None

    def __init__(self, cluster_config:dict, out_dir:PathLike = '.', prefix_name:str = '', snake_executor_file:str = None) -> None:
        self.cluster_config = cluster_config
        self.out_dir = os.path.abspath(out_dir)
        self.prefix_name = prefix_name
        if self.prefix_name: self.prefix_name+='.'
        if snake_executor_file:
            self.snake_executor_file = os.path.join(self.out_dir, snake_executor_file)
        else:
             self.snake_executor_file = snake_executor_file

        self.__cluster_validation__()

    @abstractmethod
    def __cluster_validation__(self):
        """Each scheduler should validate if the necessary options, as partition, CPUs, etc are in cluster_config.
        """

    @abstractmethod
    def build_snakemake(self):
        pass

    @abstractmethod
    def submit(self):
        pass

    def __get_full_data(self):
        data = {
            "submit_command": self.__class__.submit_command,
            "cancel_command": self.__class__.cancel_command,
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
    
    def __repr__(self):
        return f"{self.__class__.__name__}(\n{json.dumps(self.__get_full_data(), indent=5)}\n)"


class SlurmScheduler(Scheduler):
    # Override class variables
    submit_command = "sbatch"
    cancel_command = "scancel"
    shebang = "#!/bin/bash"
    job_keyword = "#SBATCH"


    def __init__(self, cluster_config: dict, out_dir: PathLike = '.', prefix_name: str = '', snake_executor_file:str = None) -> None:
        super().__init__(cluster_config = cluster_config, out_dir = out_dir, prefix_name = prefix_name, snake_executor_file = snake_executor_file)
        self.__update_internal_sbatch_values__()

    def __cluster_validation__(self):
        self.cluster_config = slurm_validation(self.cluster_config)

    def __update_internal_sbatch_values__(self):
        """This will update self.cluster_config keywords: cpus-per-task, job-name, output and error
        for better interaction with snakemake rules.
        """
        # Make log directory on demand
        cluster_log_path = os.path.join(self.out_dir, 'slurm_logs')
        cluster_log_path = os.path.abspath(cluster_log_path)
        tools.makedirs(cluster_log_path)
        # Make a copy of the user defined cluster configuration
        self._user_cluster_config = self.cluster_config.copy()
        # Update with internal values
        # threads, rule and jobid are identified and accessible during snakemake execution
        self.cluster_config.update(
            {   
                # Always use the threads defined on the rules
                "cpus-per-task": "{threads}",
                # Clear naming
                "job-name": f"{self.prefix_name}{{rule}}.{{jobid}}",
                "output": os.path.join(cluster_log_path, f"{self.prefix_name}{{rule}}.{{jobid}}.out"),
                "error": os.path.join(cluster_log_path, f"{self.prefix_name}{{rule}}.{{jobid}}.err"),
            }
        )

    def build_snakemake(self, jobs:int = 100000, latency_wait:int = 360,
                      verbose:bool = False, debug_dag:bool = False,
                      rerun_incomplete:bool = True, keep_going: bool = True) -> str:
        """Build the snakemake command

        Parameters
        ----------
        jobs : int, optional
            Use at most N CPU cluster/cloud jobs in parallel. For local execution this is an alias for --cores. Note: Set to 'unlimited' in case, this does not play a role.
            For cluster this is just a limitation.
            It is advise to provided a big number in order to do not wait for finishing of the jobs rather that launch 
            all in the queue, by default 100000
        latency_wait : int, optional
            Wait given seconds if an output file of a job is not present after the job finished. This helps if your filesystem suffers from latency, by default 120
        verbose : bool, optional
            Print debugging output, by default False
        debug_dag : bool, optional
            Print candidate and selected jobs (including their wildcards) while inferring DAG. This can help to debug unexpected DAG topology or errors, by default False
        rerun_incomplete : bool, optional
            Re-run all jobs the output of which is recognized as incomplete, by default True
        keep_going : bool, optional
            Go on with independent jobs if a job fails, by default True
        Returns
        -------
        str
            The snakemake command string.
            It also will set self._snakemake_str_cmd to the command string value
        """
        # TODO, For DEBUG Only
        if 'abfe_debug' in os.environ:
            if os.environ['abfe_debug'] == 'True':
                verbose = True
                debug_dag = True
                keep_going = False
        command = f"snakemake --jobs {jobs} --latency-wait {latency_wait} --cluster-cancel {self.cancel_command} "
        if verbose: command += "--verbose "
        if debug_dag: command += "--debug-dag "
        if rerun_incomplete: command += "--rerun-incomplete "
        if keep_going: command += "--keep-going "
        # Construct the cluster configuration
        command += f"--cluster '{self.submit_command}"
        for key in self.cluster_config:
            command += f" --{key}={self.cluster_config[key]}"
        command += "'"
        
        # Just save the command in the class
        self._snakemake_str_cmd = command

        if self.snake_executor_file:
            with open(os.path.join(self.out_dir, self.snake_executor_file), 'w') as f:
                f.write(command)
            os.chmod(os.path.join(self.out_dir, self.snake_executor_file), stat.S_IRWXU + stat.S_IRGRP + stat.S_IXGRP + stat.S_IROTH + stat.S_IXOTH)
        return command

    def submit(self, new_cluster_config:dict = None, only_build:bool = False) -> str:
        """Used to submit to the cluster the created job

        Parameters
        ----------
        new_cluster_config : dict, optional
            New definition of the cluster. It could be useful to run the snakemake command with different resources 
            as the one used on the workflow. For example, if the cluster has two partition deflt and long with 2 and 5 days as
            maximum time, we could run in the long partition the snakemake job and only ask for 1 CPU and in deflt
            the computational expensive calculations. If nothing is provided, cluster_config (passed during initialization)
            will be used, by default None
        only_build : bool, optional
            Only create the file to submit to the cluster but it will not be executed, by default False

        Returns
        -------
        str
            The output of the submit command or None.
        Raises
        ------
        RuntimeError
            If snake_executor_file is not present. You must declare it during initialization
        """
        # If extra_cluster_config, modify  self.snake_executor_file
        # Validate
        # TODO: Maybe is a good idea, instead of use the whole new_cluster_config, update the current self._user_cluster_config
        # and then validate with slurm_validation
        if new_cluster_config:
            cluster_to_work = slurm_validation(new_cluster_config)
        else:
            cluster_to_work = self._user_cluster_config
        
        # Update some configurations:
        # Make log directory on demand
        cluster_log_path = os.path.join(self.out_dir, 'slurm_logs')
        cluster_log_path = os.path.abspath(cluster_log_path)
        tools.makedirs(cluster_log_path)
        cluster_to_work.update({
            # Clear naming
            "job-name": "RuleThemAll",
            "output": os.path.join(cluster_log_path, "RuleThemAll.out"),
            "error": os.path.join(cluster_log_path, "RuleThemAll.err"),
        })

        # Create the sbatch section of the script
        sbatch_section = f"{self.shebang}\n"
        for key in cluster_to_work:
            sbatch_section += f"{self.job_keyword} --{key}={cluster_to_work[key]}\n"

        if self.snake_executor_file:
            # Update snake_executor_file
            with open(self.snake_executor_file, 'w') as sef:
                sef.write(sbatch_section + self._snakemake_str_cmd)
            if not only_build:
                # Submit to the cluster
                process = tools.run(f"{self.submit_command} {self.snake_executor_file}")
                return process.stdout
        else:
            raise RuntimeError("'snake_executor_file' attribute is not present on the current instance. Consider to call build_snakemake first")


def slurm_validation(cluster_config:dict) -> dict:
    """Validate the provided user slurm keywords

    Parameters
    ----------
    cluster_config : dict
        A dictionary with key[SBATCH keyword]: value[SBATCH value]

    Returns
    -------
    dict
        Corrected dictionary. Keywords like: c or p are translated to cpu-per-task and partition respectively.

    Raises
    ------
    ValueError
         Invalid Slurm keywords
    ValueError
        It was not provided necessary Slurm keywords
    """
    # Translate scheduler_directives
    translated_cluster_config = {}
    for key in cluster_config:
        if key not in _SBATCH_KEYWORDS: raise ValueError(f"{key} is not a valid SLURM string key")
        # Check for SBATCH flags (setting by using a boolean as value)
        if isinstance(cluster_config[key],bool):
            if cluster_config[key]:
                # Just set the flag
                translated_cluster_config[_SBATCH_KEYWORDS[key]] = ""
        else:
            translated_cluster_config[_SBATCH_KEYWORDS[key]] = cluster_config[key]
    
    # Check for important missing cluster definitions
    # TODO, check for other kwargs
    if 'partition' not in translated_cluster_config:
        raise ValueError(f"cluster_config does not have a valid SLURM definition for partition, consider to include 'p' or 'partition'")

    return translated_cluster_config

class FrontEnd:
    # TODO build a class to execute the workflow in a frontend like environment, E.g LAPTOP.
    def __init__(self) -> None:
        raise NotImplemented

def create_scheduler(scheduler_type:str, **kwargs) -> Scheduler:
    """Factory method to create the appropriate scheduler instance

    Parameters
    ----------
    scheduler_type : str
        Name of the scheduler, e.g slurm

    Returns
    -------
    Scheduler
        An instance of the proper selected scheduler

    Raises
    ------
    NotImplementedError
        In case of non implemented scheduler.
    """
    scheduler_type = scheduler_type.lower()
    if scheduler_type == "slurm":
        return SlurmScheduler(**kwargs)
    else:
        raise NotImplementedError("Invalid scheduler type. Choose from: [slurm].")

if __name__ == "__main__":
    cluster_config = {
        'p': 'deflt'
    }
    s = SlurmScheduler(cluster_config=cluster_config, out_dir='.', prefix_name='lig1', snake_executor_file='pepe.sh')
    # print(s)
    s.build_snakemake()
    s.submit()