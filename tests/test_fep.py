#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pytest


@pytest.mark.filterwarnings("ignore")
def test_fep():
    import tarfile
    import tempfile
    # import pytest
    from multiprocessing import cpu_count
    from pathlib import Path

    import yaml

    from bindflow.home import home
    from bindflow.orchestration.generate_scheduler import FrontEnd
    from bindflow.runners import calculate

    with tempfile.TemporaryDirectory(dir='.', prefix='.test_fep_') as tmp:
        home_path = home(dataDir='ci_systems')
        fname = home(dataDir='ci_systems') / 'WP6.tar.gz'
        tar = tarfile.open(fname, "r:gz")
        tar.extractall(tmp)
        tar.close()

        tmp_path = Path(tmp)/"WP6"
        ligand_files = list((tmp_path/"guest").rglob("*sdf"))[:2]

        ligands = []
        for ligand_file in ligand_files:
            ligands.append({
                'conf': ligand_file,
                'ff': {
                    'type': 'openff'
                    # 'type': 'espaloma',
                    # 'code': 'espaloma-0.3.1'
                }
            })

        protein = {
            'conf': str(tmp_path / 'host/WP6.gro'),
            'top': str(tmp_path / 'host/WP6.top'),
            'ff': {
                'code': 'espaloma-0.3.1',
            },
        }

        with open(home_path / "config-fep.yml", "r") as c:
            global_config = yaml.safe_load(c)
            # TODO
            # This is needed for MacOS when GROMACS is build wth -DGMX_GPU=OpenCL
            # This is not needed in the cluster because CUDA is different.
            global_config['extra_directives']['mdrun']['all']['ntmpi'] = 1

        num_jobs = cpu_count()
        threads = min(4, num_jobs)
        calculate(
            calculation_type='fep',
            protein=protein,
            ligands=ligands,
            membrane=None,
            cofactor=None,
            cofactor_on_protein=True,
            water_model='amber/tip3p',
            host_name='WP6',
            host_selection='resname WP6',
            hmr_factor=3,
            dt_max=0.004,
            threads=threads,
            num_jobs=num_jobs,
            replicas=1,
            scheduler_class=FrontEnd,
            job_prefix='host_guest.test',
            debug=True,
            out_root_folder_path=str(tmp_path / "fep-frontend"),
            submit=True,
            global_config=global_config)


if __name__ == '__main__':
    pass
