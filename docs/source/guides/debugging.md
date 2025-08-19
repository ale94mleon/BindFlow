(debugging-bindflow-runs)=
# Debugging BindFlow runs

To troubleshoot a BindFlow run:

1. Check the main log file
   - Located in the `.snakemake` directory of the running simulation.
   - Search for the keyword error.
   - Note the rule name and job ID associated with the failure.
2. If running on an HPC system
   - Go to the `slurm_logs` directory.
   - Find the corresponding `.err` and/or `.out` files using the rule name and job ID.
   - Review their contents — in most cases this will clearly indicate what went wrong and which steps were involved.
3. If the issue remains unclear
   - Return to the main log file and locate the reported working directory.
   - Explore that directory for more detailed logs.
   - For GROMACS simulations, check for `.log` and `.lis` files in the specified directory.
