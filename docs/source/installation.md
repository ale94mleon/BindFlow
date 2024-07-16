# 💿 Installation

## Conda dependencies

We highly recommend the latest version of [mamba](https://mamba.readthedocs.io/en/latest/index.html) for the conda environment creation. This is going to ease the resolution of dependencies. See the [Official Mamba Installation Instructions](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html). [Micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) can also be used.


`````{admonition} environment.yml
:class: tip

````{tab} Linux
  ```yaml
name: BindFlow
channels:
    - conda-forge
dependencies:
    - python >=3.8,<3.11
    - pip
    - bioconda::snakemake
    - ambertools =22.5
    - openmm =8.1.1
    - openmmforcefields =0.11.2
    - pdbfixer
    - openff-toolkit =0.14.5
    - espaloma =0.3.1
    - dglteam::dgl
  ```
````

````{tab} MacOS
  ```yaml
name: BindFlow
channels:
  - conda-forge
dependencies:
  - python >=3.8,<3.11
  - pip
    - bioconda::snakemake
    - ambertools
    - openmm
    - openmmforcefields
    - pdbfixer
    - openff-toolkit
  ```
  
  ```{note}
    In this case, users with Metal chips might need to install `espaloma =0.3.1` separately after installing [DGL](https://www.dgl.ai/pages/start.html) (probably from source). This is only necessary if the small molecule force field [Espaloma](https://github.com/choderalab/espaloma) is required. Its absence does not affect BindFlow's functionalities. GAFF will is not accessible either at the moment; see: [# 327](https://github.com/openmm/openmmforcefields/issues/327). In other words, at the moment (2024.07.12) only [OpenFF force fields](https://openforcefield.org/force-fields/force-fields/) for small molecules are accessible for MacOS users with Metal chips.
    
    This limitation is only for the automatic generation of small molecules topologies through [toff](https://toff.readthedocs.io/en/latest/index.html). User defined topologies of any kind can always be used as input for BindFlow.
  ```
````

`````

```python
  mamba env create -f environment.yml
  conda activate BindFlow
```

### Extra dependencies

```bash
conda install BindFlow
mamba install -c conda-forge "mpi4py <=3.1.5" graphviz
```

* `mpi4py <=3.1.5"`: is a dependency of [gmx_MMPBSA](https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/). used for MM(P/G)BSA calculations.
* `graphviz`: is used to create the image of the directed acyclic graph (DAG) of the workflow.

## GROMACS

[GROMACS](https://www.gromacs.org/) can be installed in various ways depending on your computer architecture. It is essential to ensure a proper installation that fits your resources. BindFlow relies on GROMACS as its molecular dynamics engine, so it is crucial to have GROMACS installed and ready to use.

Here, we will demonstrate how to build GROMACS from Source (courtesy of [Maciej Wójcik](https://biophys.uni-saarland.de/author/maciej-wojcik/)). If this method does not work, consult the [GROMACS Installation Guide](https://manual.gromacs.org/current/install-guide/index.html) for more information.

BindFlow depends on [MDAnalysis](https://www.mdanalysis.org). Current MDAnalysis-2.7.0 (2024.07.16), does not read TPR files generated for GROMACS >= 2023. So, GROMACS==2024.6 is a good option.

````{tab} Linux
  ```bash
  VERSION="2024.6"
  TARGET_LOCATION="gromacs/${VERSION}"
  SOURCE="https://gitlab.com/gromacs/gromacs.git"
  SOURCE_REF="v${VERSION}"

  mkdir -p ${TARGET_LOCATION}

  git clone --depth 1 --branch ${SOURCE_REF} "${SOURCE}" "${TARGET_LOCATION}-src"

  cmake -DGMX_GPU="CUDA" -DCMAKE_C_COMPILER=gcc-13 -DCMAKE_CXX_COMPILER=g++-13 \
      -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_INSTALL_PREFIX="$(pwd)/${TARGET_LOCATION}" -S "${TARGET_LOCATION}-src" -B "${TARGET_LOCATION}-build"

  nice -19 cmake --build "${TARGET_LOCATION}-build" --target install -j 8

  rm -rf "${TARGET_LOCATION}-build"
  rm -rf "${TARGET_LOCATION}-src"

  source "${TARGET_LOCATION}/bin/GMXRC.bash"
  ```
````

````{tab} MacOS
  Assuming [brew](https://brew.sh) is installed (a must for Mac developers).
  
  ```bash
  brew install hwloc cmake gcc@13
  ```
  
  ```bash
  VERSION=2024.6

  git clone --depth 1 --branch v${VERSION} https://gitlab.com/gromacs/gromacs.git gromacs-src
  cmake -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_INSTALL_PREFIX="$(pwd)/gromacs-${VERSION}" -DGMX_GPU=OpenCL -DGMX_HWLOC=ON -DCMAKE_C_COMPILER=gcc-13 -DCMAKE_CXX_COMPILER=g++-13 -B gromacs-build -S gromacs-src
  cmake --build "gromacs-build" --target install -j $(sysctl -n hw.logicalcpu)

  rm -rf gromacs-build
  rm -rf gromacs-src

  source gromacs-${VERSION}/bin/GMXRC.bash
  ```
  
  ```{important}
    To use the GPU, it is needed to set the following environmental variable:
    ```bash
    export GMX_GPU_DISABLE_COMPATIBILITY_CHECK=1
    ```
    Highly discouraged for production!
  ```
````

## Pip dependencies

`````{tab} With MM(P/G)BSA capabilities

  ````{tab} Production mode
    ```bash
    conda activate BindFlow
    python -m pip install -U git+https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA.git --no-deps
    git clone --depth 1 git@github.com:ale94mleon/BindFlow.git
    cd BindFlow 
    python -m pip install -e .
    ```
  ````
  ````{tab} Developer mode
    ```bash
    conda activate BindFlow
    python -m pip install -U git+https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA.git --no-deps
    git clone --depth 1 git@github.com:ale94mleon/BindFlow.git
    cd BindFlow 
    python -m pip install -e .
    ```
  ````
  
  ```{note}
    Important changes on the master branch of gmx_MMPBSA are not yet released (2024.07.16) to PyPi, so we are installing directly from the GitHub repo.
  ```
  
`````

`````{tab} Without MM(P/G)BSA capabilities

  ````{tab} Production mode
    ```bash
    conda activate BindFlow
    git clone --depth 1 git@github.com:ale94mleon/BindFlow.git
    cd BindFlow 
    python -m pip install -e .
    ```
  ````
  ````{tab} Developer mode
    ```bash
    conda activate BindFlow
    git clone --depth 1 git@github.com:ale94mleon/BindFlow.git
    cd BindFlow 
    python -m pip install -e .
    ```
  ````
`````

## Documentation dependencies

This project has an [Spinx](https://www.sphinx-doc.org/en/master/) documentation that can be built and accessed locally.

````{admonition} requirements.txt
:class: tip

  ```
myst-nb
myst-parser
sphinx_book_theme
sphinx==7.2.6
sphinx_design
sphinxcontrib-katex
sphinx-inline-tabs
sphinx_copybutton
sphinx-autobuild
  ```
````

```bash
pip install -r requirements_docs.txt
```

```bash
sphinx-autobuild docs public -a
```

Open [http://localhost:8000](http://localhost:8000). The HTML documentation is in the `public` directory.
