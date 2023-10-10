#!/usr/bin/env python3
# -*- coding: utf-8 -*-


def test_abfe():
    from abfe.home import home
    import tempfile, os, tarfile, glob, yaml
    from pathlib import Path
    # import pytest
    from multiprocessing import cpu_count
    from abfe import calculate_abfe


    with tempfile.TemporaryDirectory(dir='.', prefix='.test_abfe_') as tmp:
        home_path = Path(home(dataDir='ci_systems'))
        fname = Path(home(dataDir='ci_systems')) / 'WP6_G1.tar.gz'
        tar = tarfile.open(fname, "r:gz")
        tar.extractall(tmp)
        tar.close()
        print(os.listdir(Path(tmp) / "WP6_G1"))

        tmp_path = Path(tmp) / "WP6_G1"


        ligand_files = glob.glob(str(tmp_path / "guest/*sdf"))[:1]

        ligands = []
        for ligand_file in ligand_files:
            ligands.append({
                'conf': ligand_file,
                'ff':{
                    'type': 'espaloma',
                    'code': 'espaloma-0.3.1'
                }
            })

        protein={
            'conf': str(tmp_path / 'host/WP6.gro'),
            'top': str(tmp_path / 'host/WP6.top'),
            'ff': {
                'code': 'espaloma-0.3.1',
            },
        }

        with open(home_path / "config.yml", "r") as c:
            global_config = yaml.safe_load(c)


        # This are for the boresch restraint
        os.environ['abfe_debug_host_name'] = 'WP6'
        os.environ['abfe_debug_host_selection'] = 'resname WP6'

        num_jobs = cpu_count()
        threads = min(4, num_jobs)
        calculate_abfe(
            protein=protein,
            ligands=ligands,
            out_root_folder_path=str(tmp_path / "abfe-frontend"),
            cofactor= None,
            cofactor_on_protein=True,
            membrane = None,
            water_model = 'amber/tip3p',
            hmr_factor = 3,
            dt_max= 0.004,
            threads = threads,
            num_jobs = num_jobs,
            replicas = 1,
            submit= True,
            debug=True,
            job_prefix='host_guest.test',
            global_config = global_config)

if __name__ == '__main__':
    test_abfe()