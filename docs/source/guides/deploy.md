# BindFlow's deploy

(bindflow-deploy)=

As with any other Snakemake workflow, BindFlow can be deployed in a variety of environments, the only thing to change is how we call the [snakemake command](https://snakemake.readthedocs.io/en/stable/executing/cli.html).

At the `approach` directory, you can:

`````{admonition} Execution options
:class: tip

````{tab} FrontEnd
```bash
snakemake --jobs 12 --latency-wait 360 --rerun-incomplete --keep-incomplete --keep-going
```

All jobs will run on the current frontend. This setup is practical for testing and, in some cases, for processing a small set of ligands during MM(PB/GB)SA calculations.

To optimize the use of your frontend resources, you need to configure the `threads` (specified as a keyword to the runner) and `--jobs` (the maximum number of concurrent tasks) settings appropriately. In Snakemake terminology, `--jobs` represents the maximum number of CPUs available to run the workflow.

For example, if you set `threads = 4` and your frontend has 12 CPUs, this configuration means that for rules requiring the threads definition, you can run a maximum of 3 concurrent instances of those rules.

````

````{tab} SLURM
```bash
snakemake --jobs 100000 --latency-wait 360 --cluster-cancel scancel --rerun-incomplete --keep-incomplete --keep-going --cluster 'sbatch --partition=uds-hub --time=2-00:00:00 --gpus=1 --gres=gpu:1 --mem=4G --cpus-per-task={threads} --job-name=test.{rule}.{jobid} --output=approach/slurm_logs/test.{rule}.{jobid}.out --error=approach/slurm_logs/test.{rule}.{jobid}.err'
```

All jobs will be launched using the `sbatch` command. Note that the only differences with the frontend executions are `--cluster-cancel scancel` and `--cluster (...)`. The first defines the command used in the cluster to cancel jobs, and the second defines how to interact with the cluster. The number of jobs can be typically a big number, but there are some clusters that set a maximum number of job to be lunch, so you can lay with this number.

In this case, `--jobs` represents how many concurrent jobs can be run in parallel. This includes all jobs added to the queue, regardless of their status (RUNNING, PENDING, CANCELING, CONFIGURING, etc.). Resource allocation is handled by [SLURM](https://slurm.schedmd.com/documentation.html).

````
`````

Similar commands are generated in the `job.sh` script in the `approach` directory when the scheduler classes {py:class}`bindflow.orchestration.generate_scheduler.FrontEnd` or {py:class}`bindflow.orchestration.generate_scheduler.SlurmScheduler` is passed to any of the [BindFlow's runners](#bindflow-runners).

However, users can adapt the command based on their needs and available resources either manually (not too great 🤨) or by using the Abstract Base Class {py:class}`bindflow.orchestration.generate_scheduler.Scheduler` as a template. You can consult the source code of any of the previous schedulers to get an idea of how to implement your own.

A rough estimation of the maximum level of parallelization for the computationally intense task (those that run the GROMACS simulations) is:

* **Equilibration:** `number of ligands` X `replicas`
* **FEP:** `total number of lambda points` X `number of ligands` X `replicas`
* **MM(PB/GB)SA:** `samples` X `number of ligands` X `replicas`

These estimates only represent the parallelization at the job level. Each of these jobs is also parallelized at the GMX level, where GPU acceleration can be utilized, just like any typical GROMACS simulation (and this is even cooler 😎).

## Running workflow partially

In some cases, it may be convenient to run the workflow up to a specific point and resume it at a later time. Snakemake provides several options to achieve this.

````{tab} until

You can use the `--until` option to run the workflow up to a specific rule. This is useful when you want to stop at a particular step and complete the remaining workflow later on, potentially on different resources. For example:

```bash
snakemake (...) --until run_gmx_mmpbsa (...)
```
````

````{tab} target-jobs

Another approach is to specify target jobs by providing the rule name and the associated wildcards. This allows you to execute the dependencies of a specific job and its instances. For example:

```bash
snakemake (...) --target-jobs run_gmx_mmpbsa:ligand_name=ligand1,replica=3,sample=2 (...)
```

The above command will execute all dependent rules for `run_gmx_mmpbsa` and specifically the instances of `run_gmx_mmpbsa` with the following wildcards:

```yml
ligand_name: ligand1
replica: 3
sample: 2
```
````
