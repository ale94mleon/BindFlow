# Handling cofactors properly during MMGBSA

Modify the command from the `job.sh` adding `--until build_ligand_system`.

```bash
snakemake --until build_ligand_system --jobs 100000 --latency-wait 360 --cluster-cancel scancel --rerun-incomplete --keep-incomplete --keep-going --cluster 'sbatch --partition=uds-hub --time=0-20:00:00 --gpus=1 --gres=gpu:1 --mem=5G --ntasks={threads} --cpus-per-task=1 --job-name=a2a.espaloma-0.3.1.{rule}.{jobid} --output=/scratch/uds_alma015/smaug/data/users/alejandro/simulation/BindFlow_simulations/A2A/mmpbsa/espaloma-0.3.1/slurm_logs/a2a.espaloma-0.3.1.{rule}.{jobid}.out --error=/scratch/uds_alma015/smaug/data/users/alejandro/simulation/BindFlow_simulations/A2A/mmpbsa/espaloma-0.3.1/slurm_logs/a2a.espaloma-0.3.1.{rule}.{jobid}.err'
```

After conclusion, we used the script `make_ndx.py` to modify the RECEPTOR group such as include the cofactors (in our case Na ion for A2A).

Then, we restore the command by removing `--until build_ligand_system`, and execute.

```bash
snakemake --jobs 100000 --latency-wait 360 --cluster-cancel scancel --rerun-incomplete --keep-incomplete --keep-going --cluster 'sbatch --partition=uds-hub --time=0-20:00:00 --gpus=1 --gres=gpu:1 --mem=5G --ntasks={threads} --cpus-per-task=1 --job-name=a2a.espaloma-0.3.1.{rule}.{jobid} --output=/scratch/uds_alma015/smaug/data/users/alejandro/simulation/BindFlow_simulations/A2A/mmpbsa/espaloma-0.3.1/slurm_logs/a2a.espaloma-0.3.1.{rule}.{jobid}.out --error=/scratch/uds_alma015/smaug/data/users/alejandro/simulation/BindFlow_simulations/A2A/mmpbsa/espaloma-0.3.1/slurm_logs/a2a.espaloma-0.3.1.{rule}.{jobid}.err'
```

## Root of the problem

```python
host_name='Protein or resname yCA or resname COF'
```

Will make crash the gmx_MMPBSA because actually `center_xtc` is failing before when dealing with this non standard `host_name`.
This is related with issue [#32](https://github.com/ale94mleon/BindFlow/issues/32).
