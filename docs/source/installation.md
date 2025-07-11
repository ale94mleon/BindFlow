# 💿 Installation

## Conda dependencies

We highly recommend the latest version of [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) for the conda environment creation and a fresh environment (as it will be demonstrated here). This is going to ease the resolution of dependencies. See the [Official Micromamba Installation Instructions](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html). [Mamba](https://mamba.readthedocs.io/en/latest/index.html) may also be used.

`````{admonition} environment.yml
:class: tip

We would like to work in a more relaxed environment, but we have encountered challenging situations. For now, we recommend the following pinned environment.

````{tab} Linux 🐧

```{eval-rst}
.. literalinclude:: env/latest/environment-linux.yml
   :language: yaml
```

````

````{tab} MacOS 🍏

```{eval-rst}
.. literalinclude:: env/latest/environment-macos.yml
   :language: yaml
```

````

```{note}
We observed that conda must be configured with `channel_priority: flexible` instead of `strict`.
```


`````

```python
micromamba env create -f environment.yml --channel-priority flexible -y
micromamba activate BindFlow
```

## GROMACS

[GROMACS](https://www.gromacs.org/) can be installed in various ways depending on your computer architecture. It is essential to ensure a proper installation that fits your resources. BindFlow relies on GROMACS as its molecular dynamics engine, so it is crucial to have GROMACS installed and ready to use.

Here, we will demonstrate how to build GROMACS from Source (courtesy of [Maciej Wójcik](https://biophys.uni-saarland.de/author/maciej-wojcik/)). If this method does not work, consult the [GROMACS Installation Guide](https://manual.gromacs.org/current/install-guide/index.html) for more information.

BindFlow depends on [MDAnalysis](https://www.mdanalysis.org). Current MDAnalysis-2.7.0 (2024.07.16), does not read TPR files generated for GROMACS >= 2023. So, GROMACS==2022.6 is a good option.

````{tab} Linux 🐧
  ```bash
  VERSION="2022.6"
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

````{tab} MacOS 🍏
  Assuming [brew](https://brew.sh) is installed (a must for Mac developers).
  
  ```bash
  brew install hwloc cmake gcc@13
  ```
  
  ```bash
  VERSION=2022.4

  git clone --depth 1 --branch v${VERSION} https://gitlab.com/gromacs/gromacs.git gromacs-src
  cmake -DGMX_BUILD_OWN_FFTW=ON -DCMAKE_INSTALL_PREFIX="$(pwd)/gromacs-${VERSION}" -DGMX_GPU=OpenCL -DGMX_HWLOC=ON -DCMAKE_C_COMPILER=gcc-13 -DCMAKE_CXX_COMPILER=g++-13 -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -B gromacs-build -S gromacs-src
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
    micromamba activate BindFlow
    git clone --depth 1 git@github.com:ale94mleon/BindFlow.git
    cd BindFlow 
    python -m pip install -e . --no-deps
    python -m pip install -U git+https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA.git@27929e02067bc2321286809818d778a77a872010 --no-deps
    ```
  ````
  ````{tab} Developer mode
    ```bash
    micromamba activate BindFlow
    git clone --depth 1 git@github.com:ale94mleon/BindFlow.git
    cd BindFlow 
    python -m pip install -e . --no-deps
    python -m pip install -U git+https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA.git@27929e02067bc2321286809818d778a77a872010 --no-deps
    ```
  ````
  
  ```{note}
  Currently, we are using the `gmx_MMPBSA` commit [27929e0](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/commit/27929e02067bc2321286809818d778a77a872010). This commit has been tested and provides flexibility in selecting the Python version.
  ```
  
`````

`````{tab} Without MM(P/G)BSA capabilities

  ````{tab} Production mode
    ```bash
    micromamba activate BindFlow
    git clone --depth 1 git@github.com:ale94mleon/BindFlow.git
    cd BindFlow 
    python -m pip install -e . --no-deps
    ```
  ````
  ````{tab} Developer mode
    ```bash
    micromamba activate BindFlow
    git clone --depth 1 git@github.com:ale94mleon/BindFlow.git
    cd BindFlow 
    python -m pip install -e . --no-deps
    ```
  ````
`````

## Documentation dependencies

This project has an [Sphinx](https://www.sphinx-doc.org/en/master/) documentation that can be built and accessed locally.

````{admonition} requirements_doc.txt
:class: tip

  ```
  myst-nb
  myst-parser
  sphinx_book_theme
  sphinx==7.2.6
  sphinx_design
  sphinxcontrib-katex
  sphinxcontrib-mermaid
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

## Testing installation

Finally, it is advised to check if everything is alright. Be patient and go for a coffee ☕, this could take a couple of minutes (~11 min on my Laptop) ⏳.

````{tab} BindFlow is already cloned
  ```bash
  cd BindFlow # The path to your local copy of the repository
  python -m pytest test
  ```
````

````{tab} BindFlow is not cloned yet
  ```bash
  git clone --depth 1 git@github.com:ale94mleon/BindFlow.git
  cd BindFlow
  python -m pytest test
  ```
````

`````{admonition} Expected results
:class: info
````{tab} Linux 🐧
  ```yaml
  passed: 5
  xfailed: 0
  ```
````
````{tab} MacOS 🍏
  ```bash
  passed: 3
  xfailed: 2  # from test/test_small.py
  ```
````
`````
