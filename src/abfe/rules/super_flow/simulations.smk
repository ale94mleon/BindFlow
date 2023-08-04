# Consider to use 
from abfe import root_path as abfe_root_path
approach_path = config["out_approach_path"]


# Do Check all results
rule abfe_ligand_result:
    input:
        dG_path = expand(approach_path + "/{ligand_name}/{replica}/dG_results.csv", ligand_name = config['ligand_names'], replica = list(map(str, range(1,1 + config['replicas']))))

# Do Equilibration
include: abfe_root_path + '/rules/equi/Snakefile'

# Do FEP
include: abfe_root_path + '/rules/fep/Snakefile'