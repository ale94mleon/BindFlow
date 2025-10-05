# 🚀 Tutorials

Get started quickly with **ready-to-simulate** examples provided in the [examples](https://github.com/ale94mleon/BindFlow/tree/main/examples) directory of the BindFlow repository. These examples are designed to work out of the box, giving you a hands-on introduction to BindFlow.

## Current examples

In the example directory there are three systems:

- **A2A**: membrane protein with a Na ion cofactor,
- **CyclophilinD**: soluble protein,
- **SAMPL6-OA**: non-protein receptor, ideal for testing as is the smallest system.

For each system you will find:

- Protein and ligands input files in the `inputs` directory,
- BindFlow configuration files for:
  - FEP: `config-fep.yml`, and
  - MM(PB/GB)SA: `config-mmpbsa.yml`.
- BindFlow executor files: `executor-fep.py` and `executor-mmpbsa.py`.

Files presented on the examples are those used on the [BindFlow's paper](https://www.biorxiv.org/content/10.1101/2025.09.25.678545v1)

Here we will only follow step-by-step the **CyclophilinD** example, the other examples are very similar.

## Tips before starting

1. **Review the examples** - Browse through the example simulations to understand how they are structured.
2. **Adjust configurations** - Update the configuration files to match your cluster resources and environment. This ensures simulations run smoothly.
3. **Reference configuration files** - Several files are provided for illustration in the [config](https://github.com/ale94mleon/BindFlow/tree/main/examples/configs) directory. Use them as a guide to create or customize your own simulation setups.
4. **Run your simulations** - Once configured, execute the examples to see BindFlow in action. For deployment instructions, see the [Deploy guide](./guides/deploy.md).
5. **Troubleshooting** - Check common issues and solutions in the [FAQ](./faq.md). If your problem is not listed, questions on how to use BindFlow, or if you want to give feedback or share ideas and new features, please head to the [BindFlow's discussions](https://github.com/ale94mleon/BindFlow/discussions).

## CylophilinD FEP tutorial

### Hardware

To complete this tutorial we need an HPC with SLURM as task manager. This is not strictly needed for BindFlow to run, but it will be easy to be on the same page for the tutorial.

```{admonition} Hardware accessibility
:class: info

In case you do not have access to an HPC-like computer environment and you are only interesting in see if BindFlow works, consider to use the small SAMPL6-OA system and only a couple of ligands. So the calculations are doable in your workstation.
```

### Getting ready

First we need to install BindFlow in our computer system, for that follow the [Installation instructions](./installation.md).

Get the files by cloning the [BindFlow repo](https://github.com/ale94mleon/BindFlow/)

```bash
git clone --depth 1 git@github.com:ale94mleon/BindFlow.git
cd BindFlow/examples/CyclophilinD
```

#### Exploring and modifying the executor file

Explore the directory, check the files and try to understand the logic of `executor-fep.py`. You will see that only a single function is used: {py:func}`bindflow.runners.calculate`. This is the core function you should get familiar for most of the general cases. Therefore, it is **highly recommended** to learn its parameters and options. The documentation of {py:func}`bindflow.runners.calculate` is very extensive, but you will learn how to use BindFlow effectively.

Let's modify `executor-fep.py` by:

- using only one of the three force fields,
- using two replicas,
- using only two ligands and,
- submitting the main job to the queue system.

```diff
10,17d9
<     "openff_unconstrained-2.0.0": {
<         "type": "openff",
<         "code": "openff_unconstrained-2.0.0.offxml",
<     },
<     "espaloma-0.3.1": {
<         "type": "espaloma",
<         "code": "espaloma-0.3.1",
<     },
43c35
<         ligands=ligands,
---
>         ligands=ligands[:2],
49c41
<         replicas=3,
---
>         replicas=2,
53c45
<         submit=False,
---
>         submit=True,
```

#### Tuning the config file

The provided `config-fep.yml` offers SLURM configurations for the calculation jobs (`cluster:options:calculation`) and for the sentinel job (`cluster:options:job`). The last will only launch jobs and handled the snakemake queue system, but it will not perform heavy calculations, it will be mainly sleeping and waiting. Keywords for these sections are any valid SLURM parameter. You will need to adjust based on your HPC configuration.

`extra_directives/dependencies` is a list of commands that will be executed sequentially before any GROMACS command. In the provided example, this section is used to load GROMACS package on our system (for sure different in your case) and avoid GROMACS backup on files (this ensures cleaner working directory during production). You **must** change those sections accordingly to your HPC configuration.

#### Run BindFlow run

After you tuned your configurations, you can just:

```bash
python executor-fep.py
```

#### Checking the run

BindFlow offers some basic, but powerful CLI functionalities; check the command `bindflow -h`.

You can also check the SLURM queue and the main snakemake log file located at `fep/gaff-2.11/.snakemake/log` to see the completion percent of the pipeline.

To debug possible issues, read the [Debugging page](./guides/debugging.md) and the [FAQ](./faq.md).

```{admonition} Friendly reminder
:class: info

The examples are a great way to get hands-on quickly, but they’re **not a replacement for the full documentation**. For a deeper understanding of BindFlow and its features, we encourage you to explore the documentation--we put a lot of love into it! 😊
```
