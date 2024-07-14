Installation
============

We will use `mamba <https://mamba.readthedocs.io/en/latest/>`__ and `environment.yml <https://github.com/ale94mleon/BindFlow/blob/main/environment.yml>`__ definition.

If you do not have ``mamba`` installed, then:

.. code-block:: bash

  conda install mamba -n base -c conda-forge

.. warning::

  You could try also with `conda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/windows.html>`__ but it could take a while to build the environment. We observed that mamba builds faster

.. code-block:: bash

  mamba env create -f environment.yml

For MacOS users with ARM ships you might need to remove from the environment file:

.. code-block:: yaml
  
  - espaloma >=0.3.1
  - dglteam::dgl

``espaloma``_ requires ``dgl 1.1.2.*``, which is not for osx-arm64

If you want to modify the code and contribute, then:

.. code-block:: bash

    git clone https://github.com/ale94mleon/BindFlow.git
    cd BindFlow 
    conda activate BindFlow
    python -m pip install -e .


In order to MDAnalysis-2.7.0 read the GROMACS TPR file you must install GROMACS<2023

To get gmx_MMPBSA (the current PyPi version is not yet updated):

.. code-block:: bash
  
  python -m pip install -U git+https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA.git --no-deps 
  
