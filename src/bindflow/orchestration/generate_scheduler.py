import json
import os
import stat
from abc import ABC, abstractmethod
from pathlib import Path

from bindflow.utils import tools
from bindflow.utils.cluster import _SBATCH_KEYWORDS
from bindflow.utils.tools import PathLike


class Scheduler(ABC):
    """Abstract Base Class to build an Schedular

    Class variables
    ---------------
    submit_command : str
        The command used for your scheduler to launch jobs
    cancel_command : str
        Command used to cancel jobs
    shebang : str
        Used to build script and detect properly the environment E.g: ``#!/bin/bash``, ``#!/bin/sh``, ...
        This will be used to make the ``snake_executor_file`` executable.
    """
    # Default class variables
    submit_command = None
    cancel_command = None
    shebang = None
    job_keyword = None

    def __init__(self, cluster_config: dict, out_dir: PathLike = '.', prefix_name: str = '', snake_executor_file: str = None) -> None:
        """Constructor of the class

        Parameters
        ----------
        cluster_config : dict
            All the necessary information for the specific schedular
        out_dir : PathLike, optional
            Where all files will be exported and executed, by default '.'
        prefix_name : str, optional
            A prefix append to the jobs names for easy identification, by default ''
        snake_executor_file : str, optional
            The name/path of the file that will be used for execution of the workflow, by default None
        """
        self.cluster_config = cluster_config
        self.out_dir = Path(out_dir).resolve()
        self.prefix_name = prefix_name
        if self.prefix_name:
            self.prefix_name += '.'
        if snake_executor_file:
            self.snake_executor_file = self.out_dir/snake_executor_file
        else:
            self.snake_executor_file = snake_executor_file

        self.__cluster_validation__()

    @abstractmethod
    def __cluster_validation__(self):
        """Each scheduler should validate if the necessary options, as partition, CPUs, etc are in cluster_config.
        """
        pass

    @abstractmethod
    def build_snakemake(self, jobs: int):
        """Function to create the snakemake command

        Parameters
        ----------
        jobs : int
            Number of snakemake jobs. Passed to the flag `--jobs`
        """
        pass

    @abstractmethod
    def submit(self, new_cluster_config: dict, only_build: bool, job_prefix: str):
        """Command to submit the jobs

        Parameters
        ----------
        new_cluster_config : dict
            If an specific cluster configuration is wanted for the
            man Snakemake job (this is mainly hanging and waiting
            for rules completion)
        only_build : bool
            Only build the files but do not execute the command
        job_prefix : str
            A job prefix identification for the cluster
        """

    def __get_full_data(self) -> dict:
        """Get the data of the class

        Returns
        -------
        dict
            Information of the class
        """
        data = {
            "submit_command": self.__class__.submit_command,
            "cancel_command": self.__class__.cancel_command,
            "shebang": self.__class__.shebang,
            "job_keyword": self.__class__.job_keyword,
        }
        data.update(self.__dict__)
        return data

    def to_json(self, out_file: str = "cluster.json"):
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

    def __init__(self, cluster_config: dict, out_dir: PathLike = '.', prefix_name: str = '', snake_executor_file: str = None) -> None:
        super().__init__(cluster_config=cluster_config, out_dir=out_dir, prefix_name=prefix_name, snake_executor_file=snake_executor_file)
        self.__update_internal_sbatch_values__()

    def __cluster_validation__(self):
        self.cluster_config = slurm_validation(self.cluster_config)

    def __update_internal_sbatch_values__(self):
        """This will update self.cluster_config keywords: ntasks, cpus-per-task, job-name, output and error
        for better interaction with snakemake rules.
        """
        # Make log directory on demand
        cluster_log_path = (self.out_dir/'slurm_logs').resolve()
        cluster_log_path.mkdir(exist_ok=True, parents=True)
        # Make a copy of the user defined cluster configuration
        self._user_cluster_config = self.cluster_config.copy()
        # Update with internal values
        # threads, rule and jobid are identified and accessible during snakemake execution
        self.cluster_config.update(
            {
                # Always use the threads defined on the rules
                # Need to define in this way so MPI process detect slots properly.
                "ntasks": "{threads}",
                "cpus-per-task": "1",
                # Clear naming
                "job-name": f"{self.prefix_name}{{rule}}.{{jobid}}",
                "output": cluster_log_path/f"{self.prefix_name}{{rule}}.{{jobid}}.out",
                "error": cluster_log_path/f"{self.prefix_name}{{rule}}.{{jobid}}.err",
            }
        )

    def build_snakemake(self, jobs: int = 100000, latency_wait: int = 360,
                        verbose: bool = False, debug_dag: bool = False,
                        rerun_incomplete: bool = True, keep_incomplete: bool = True,
                        keep_going: bool = True) -> str:
        """Build the snakemake command
        TODO Consider to put it in the parent class

        Parameters
        ----------
        jobs : int, optional
            Use at most N CPU cluster/cloud jobs in parallel. For local execution this is an alias for --cores.
            Note: Set to 'unlimited' in case, this does not play a role.
            For cluster this is just a limitation.
            It is advise to provided a big number in order to do not wait for finishing of the jobs rather that launch
            all in the queue, by default 100000
        latency_wait : int, optional
            Wait given seconds if an output file of a job is not present after the job finished.
            This helps if your filesystem suffers from latency, by default 120
        verbose : bool, optional
            Print debugging output, by default False
        debug_dag : bool, optional
            Print candidate and selected jobs (including their wildcards) while inferring DAG.
            This can help to debug unexpected DAG topology or errors, by default False
        rerun_incomplete : bool, optional
            Re-run all jobs the output of which is recognized as incomplete, by default True
        keep_incomplete : bool, optional
            TODO !!! This could let to undesired effects but it is needed for GROMACS continuation
            Do not remove incomplete output files by failed jobs, by default True.
        keep_going : bool, optional
            Go on with independent jobs if a job fails, by default True
        Returns
        -------
        str
            The snakemake command string.
            It also will set self._snakemake_str_cmd to the command string value
        """
        # TODO, For DEBUG Only
        if 'BINDFLOW_DEBUG' in os.environ:
            if os.environ['BINDFLOW_DEBUG'] == 'True':
                verbose = True
                debug_dag = True
                keep_going = False
        command = f"snakemake --jobs {jobs} --latency-wait {latency_wait} --cluster-cancel {self.cancel_command} "
        if verbose:
            command += "--verbose "
        if debug_dag:
            command += "--debug-dag "
        if rerun_incomplete:
            command += "--rerun-incomplete "
        if keep_incomplete:
            command += "--keep-incomplete "
        if keep_going:
            command += "--keep-going "
        # Construct the cluster configuration
        command += f"--cluster '{self.submit_command}"
        # Here is the only possible difference, maybe it could be creates an
        # abstract method that return cluster_config to a string representation valid to execute the jobs
        for key in self.cluster_config:
            command += f" --{key}={self.cluster_config[key]}"
        command += "'"

        # Just save the command in the class
        self._snakemake_str_cmd = command

        if self.snake_executor_file:
            with open(self.out_dir/self.snake_executor_file, 'w') as f:
                f.write(command)
            os.chmod(self.out_dir/self.snake_executor_file, stat.S_IRWXU + stat.S_IRGRP + stat.S_IXGRP + stat.S_IROTH + stat.S_IXOTH)
        return command

    def submit(self, new_cluster_config: dict = None, only_build: bool = False, job_prefix: str = "") -> str:
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
        job_prefix : bool, optional
            It will be added as {job_prefix}.RuleThemAll , by default False
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
        cluster_log_path = (self.out_dir/'slurm_logs').resolve()
        cluster_log_path.mkdir(exist_ok=True, parents=True)
        cluster_to_work.update({
            # Clear naming
            "job-name": f"{job_prefix}.RuleThemAll",
            "output": cluster_log_path/f"{job_prefix}.RuleThemAll.out",
            "error": cluster_log_path/f"{job_prefix}.RuleThemAll.err",
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


def slurm_validation(cluster_config: dict) -> dict:
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
        if key not in _SBATCH_KEYWORDS:
            raise ValueError(f"{key} is not a valid SLURM string key")
        # Check for SBATCH flags (setting by using a boolean as value)
        if isinstance(cluster_config[key], bool):
            if cluster_config[key]:
                # Just set the flag
                translated_cluster_config[_SBATCH_KEYWORDS[key]] = ""
        else:
            translated_cluster_config[_SBATCH_KEYWORDS[key]] = cluster_config[key]

    # Check for important missing cluster definitions
    # TODO, check for other kwargs
    if 'partition' not in translated_cluster_config:
        raise ValueError("cluster_config does not have a valid SLURM definition for partition, consider to include 'p' or 'partition'")

    return translated_cluster_config


class FrontEnd(Scheduler):
    # Override class variables
    submit_command = "bash"
    shebang = "#!/bin/bash"

    # TODO build a class to execute the workflow in a frontend-like environment, E.g LAPTOP.
    def __init__(self, cluster_config: None = None, out_dir: PathLike = '.', prefix_name: str = '', snake_executor_file: str = None) -> None:
        super().__init__(cluster_config=cluster_config, out_dir=out_dir, prefix_name=prefix_name, snake_executor_file=snake_executor_file)

    def __cluster_validation__(self): ...

    def build_snakemake(self, jobs: int = 12, latency_wait: int = 360,
                        verbose: bool = False, debug_dag: bool = False,
                        rerun_incomplete: bool = True, keep_incomplete: bool = True,
                        keep_going: bool = True) -> str:
        """Build the snakemake command
        TODO Consider to put it in the parent class

        Parameters
        ----------
        jobs : int, optional
            Use at most N CPU cluster/cloud jobs in parallel. For local execution this is an alias for --cores.
            Note: Set to 'unlimited' in case, this does not play a role.
            For cluster this is just a limitation.
            It is advise to provided a big number in order to do not wait for finishing of the jobs rather that launch
            all in the queue, by default 100000
        latency_wait : int, optional
            Wait given seconds if an output file of a job is not present after the job finished.
            This helps if your filesystem suffers from latency, by default 120
        verbose : bool, optional
            Print debugging output, by default False
        debug_dag : bool, optional
            Print candidate and selected jobs (including their wildcards) while inferring DAG.
            This can help to debug unexpected DAG topology or errors, by default False
        rerun_incomplete : bool, optional
            Re-run all jobs the output of which is recognized as incomplete, by default True
        keep_incomplete : bool, optional
            TODO !!! This could let to undesired effects but it is needed for GROMACS continuation
            Do not remove incomplete output files by failed jobs, by default True.
        keep_going : bool, optional
            Go on with independent jobs if a job fails, by default True
        Returns
        -------
        str
            The snakemake command string.
            It also will set self._snakemake_str_cmd to the command string value
        """
        # TODO, For DEBUG Only
        if 'BINDFLOW_DEBUG' in os.environ:
            if os.environ['BINDFLOW_DEBUG'] == 'True':
                verbose = True
                debug_dag = True
                keep_going = False
        command = f"snakemake --jobs {jobs} --latency-wait {latency_wait} "
        if verbose:
            command += "--verbose "
        if debug_dag:
            command += "--debug-dag "
        if rerun_incomplete:
            command += "--rerun-incomplete "
        if keep_incomplete:
            command += "--keep-incomplete "
        if keep_going:
            command += "--keep-going "

        # Just save the command in the class
        self._snakemake_str_cmd = command

        if self.snake_executor_file:
            with open(self.out_dir/self.snake_executor_file, 'w') as f:
                f.write(command)
            os.chmod(self.out_dir/self.snake_executor_file, stat.S_IRWXU + stat.S_IRGRP + stat.S_IXGRP + stat.S_IROTH + stat.S_IXOTH)
        return command

    def submit(self, only_build: bool = False, **kwargs) -> str:
        """Used to submit to the cluster the created job

        Parameters
        ----------
        only_build : bool, optional
            Only create the file to submit to the Frontend but it will not be executed, by default False
        **kwargs : object, optional
            This is only added for compatibility. t is like this in order that the snakemake rules can pass arguments irrespective if is
            SLURM FrontEnd without to check for configuration.
            In reality it will not be used at all on this method.
        Returns
        -------
        str
            The output of the submit command or None.
        Raises
        ------
        RuntimeError
            If snake_executor_file is not present. You must declare it during initialization
        """
        # Create the sbatch section of the script
        bash_section = f"{self.shebang}\n"

        if self.snake_executor_file:
            # Update snake_executor_file
            with open(self.snake_executor_file, 'w') as sef:
                sef.write(bash_section + self._snakemake_str_cmd)
            if not only_build:
                # Submit to the Frontend
                tools.run(f"{self.submit_command} {self.snake_executor_file}", interactive=True)
        else:
            raise RuntimeError("'snake_executor_file' attribute is not present on the current instance. Consider to call build_snakemake first")


if __name__ == "__main__":
    pass
