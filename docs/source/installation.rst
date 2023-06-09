Installation
============

We will use `mamba <https://mamba.readthedocs.io/en/latest/>`__. First, you must download `environment.yml <https://github.com/ale94mleon/BindFlow/blob/main/environment.yml>`__.
In case some problems happens, you can use `environment.yml <https://github.com/ale94mleon/BindFlow/blob/main/environment_pinned.yml>`__

If you do not have ``mamba`` installed, then:

.. code-block:: bash

  conda install mamba -n base -c conda-forge

.. warning::

  You could try also with `conda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/windows.html>`__ but it could take a while building the environment.
  We observed that mamba build faster the environment; however, it could also take some time solving the dependencies. Please, be patient. If some error happens, then use the
  `environment_pinned.yml <https://github.com/bigginlab/ABFE_workflow/blob/main/environment_pinned.yml>`__ instead.

.. code-block:: bash

  mamba env create -f environment.yml

If you want to modify the code and contribute, then:

.. code-block:: bash

    git clone https://github.com/bigginlab/ABFE_workflow.git
    cd ABFE_workflow 
    conda activate abfe
    python -m pip install -e .
