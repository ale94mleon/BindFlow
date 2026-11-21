# 💿 Installation

BindFlow is installed in **two** conda environments, created and linked by the BindFlow installer:

- **`BindFlow`**: the workflow, the free energy analysis and (optionally) MM(PB/GB)SA. This is the environment you activate and work with.
- **`BindFlow-TOFF`**: small-molecule parameterization ([TOFF](https://github.com/ale94mleon/TOFF): OpenFF, GAFF and optionally Espaloma) and [pdbfixer](https://github.com/openmm/pdbfixer). You never need to activate it; BindFlow calls it when needed.

Keeping the heavy OpenMM/OpenFF/PyTorch stack in its own environment makes the dependency resolution much faster and more robust.

```{warning}
**macOS Intel (x86_64) is not supported** and the installer stops there. On that platform, the conda-forge build of AmberTools' `antechamber` randomly reads zero atoms from valid SDF files, so small-molecule parameterization with OpenFF (AM1-BCC charges) and GAFF fails intermittently. This is an upstream bug (see the *Known issues* in the [TOFF changelog](https://github.com/ale94mleon/TOFF/blob/main/docs/source/CHANGELOG.md)). Supported platforms are Linux (x86_64 and ARM) and macOS Apple Silicon. On Apple Silicon, run the installer from a native shell, not under Rosetta.
```

## Installing micromamba

The installer needs [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) (recommended) or [mamba](https://mamba.readthedocs.io/en/latest/index.html) >= 2. See the [Official Micromamba Installation Instructions](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html).

## Running the installer

Download the installer and choose the options you need:

```bash
curl -fsSL -O https://raw.githubusercontent.com/ale94mleon/BindFlow/main/install/installer.sh
```

`````{tab} Without MM(P/G)BSA capabilities

  ````{tab} 🟢 Stable Release
    ```bash
    bash installer.sh
    ```
  ````
  ````{tab} 🟡 Install from GitHub (latest commit)
    ```bash
    bash installer.sh --bindflow git+https://github.com/ale94mleon/BindFlow.git
    ```
  ````
  ````{tab} 🔵 Developer Mode (editable)
    ```bash
    git clone --depth 1 git@github.com:ale94mleon/BindFlow.git
    bash BindFlow/install/installer.sh --editable BindFlow
    ```
  ````
`````

`````{tab} With MM(P/G)BSA capabilities

  ````{tab} 🟢 Stable Release
    ```bash
    bash installer.sh --mmpbsa
    ```
  ````
  ````{tab} 🟡 Install from GitHub (latest commit)
    ```bash
    bash installer.sh --mmpbsa --bindflow git+https://github.com/ale94mleon/BindFlow.git
    ```
  ````
  ````{tab} 🔵 Developer Mode (editable)
    ```bash
    git clone --depth 1 git@github.com:ale94mleon/BindFlow.git
    bash BindFlow/install/installer.sh --mmpbsa --editable BindFlow
    ```
  ````

  ```{note}
  Currently, we are using the `gmx_MMPBSA` commit [27929e0](https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA/commit/27929e02067bc2321286809818d778a77a872010). This commit has been tested and provides flexibility in selecting the Python version.
  ```
`````

Add `--espaloma` to any of the commands above to also install the [Espaloma](https://github.com/choderalab/espaloma) force field.

```{important}
Espaloma depends on [DGL](https://www.dgl.ai), which conda-forge does not provide for **ARM Linux** (`linux-aarch64`, e.g. Grace/Graviton nodes or Docker on Apple Silicon). It is available on x86_64 Linux and on macOS Apple Silicon.
```

The installer ends with a quick check (it parameterizes ethane with every installed force field). When it finishes, activate BindFlow with:

```bash
micromamba activate BindFlow
```

````{admonition} All installer options
:class: tip

  ```text
  --mmpbsa           Add MM(PB/GB)SA support (AmberTools + gmx_MMPBSA) to the main environment
  --espaloma         Add Espaloma to the parameterization environment
  --name NAME        Main environment name (default: BindFlow)
  --toff-name NAME   Parameterization environment name (default: BindFlow-TOFF)
  --bindflow SPEC    pip spec for BindFlow (default: bindflow)
  --editable DIR     Install BindFlow in editable mode from a local clone instead
  --toff SPEC        pip spec for TOFF (default: toff>=0.3.0)
  --force            Remove existing environments with the same names first
  -h, --help         Show this help
  ```
````

### How BindFlow finds `BindFlow-TOFF`

BindFlow looks for the parameterization environment in this order:

1. The `toff_prefix` argument of {py:func}`bindflow.runners.calculate`.
2. The `BINDFLOW_TOFF_PREFIX` environment variable. The installer sets it every time `BindFlow` is activated.
3. An environment called `BindFlow-TOFF` next to the `BindFlow` environment.

If you used the default names, you do not need to do anything. If you are running on a cluster, make sure `BindFlow-TOFF` is on a file system that the compute nodes can also access.

## GROMACS

[GROMACS](https://www.gromacs.org/) can be installed in various ways depending on your computer architecture. It is essential to ensure a proper installation that fits your resources. BindFlow relies on GROMACS as its molecular dynamics engine, so it is crucial to have GROMACS installed and ready to use.

Here, we will demonstrate how to build GROMACS from Source (courtesy of [Maciej Wójcik](https://biophys.uni-saarland.de/author/maciej-wojcik/)). If this method does not work, consult the [GROMACS Installation Guide](https://manual.gromacs.org/current/install-guide/index.html) for more information.

```{important}
At the moment we support GROMACS versions in the range [2022, 2027); BindFlow throws a clean error at the very beginning of the execution of {py:func}`bindflow.runners.calculate` if otherwise.
```

<!--
2025 -> TPX 137 (until MDAnalysis 2.10)
2026 -> TPX 138
 
Check the TPX version of the new GROMACS:

1. Create a TPR with the new GROMACS
2. Then execute:

```python
import MDAnalysis as mda
from MDAnalysis.topology import tpr
from MDAnalysis import lib

with lib.util.openany("prod.tpr", mode="rb") as infile:
    tprf = infile.read()
data = tpr.utils.TPXUnpacker(tprf)
tpr.utils.read_tpxheader(data)
```
IF the following pass, is ok. You can aslo, instead the last command:

```python
data.do_string()
data.unpack_int()
print(data.unpack_int(), "this is the TPR")
```
And check if the version is supported by https://docs.mdanalysis.org/stable/documentation_pages/topology/tpr_util.html#module-MDAnalysis.topology.tpr. Maybe a new version of MDAnalysis.-->

<!-- ```{important}
BindFlow depends on [MDAnalysis](https://www.mdanalysis.org). Current `MDAnalysis-2.9.0` (2026.03.31), does not read TPR files generated by `GROMACS >= 2026`; so a lower version is needed. `GROMACS==2025.4` is a good option.

If `GROMACS >= 2023`, BindFlow throws a clean error at the very beginning of the execution of {py:func}`bindflow.runners.calculate`.
``` -->

````{tab} Linux 🐧
  ```bash
  VERSION="2025.4"
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
  VERSION=2025.4

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

## Documentation dependencies

This online [Sphinx](https://www.sphinx-doc.org/en/master/) documentation can also be built and accessed locally.

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
  roman
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

```bash
micromamba activate BindFlow
micromamba install -c conda-forge pytest
pytest --pyargs bindflow.tests
```

```{admonition} Expected results
:class: info

- `test_fep`, `test_mmpbsa` and the OpenFF parameterization test should pass. `test_mmpbsa` needs the `--mmpbsa` option of the installer.
- The Espaloma test is skipped if Espaloma is not installed.
- On macOS, the GAFF and Espaloma tests are marked as expected to fail (`xfail`).
```

```{note}
The installation is tested automatically on GitHub Actions on Linux (x86_64 and ARM) and macOS Apple Silicon on every change of the installer or the code. To keep this affordable, those runs only include the installer check and the small-molecule parameterization tests; `test_fep` and `test_mmpbsa` run real GROMACS simulations and are only run locally, before each release.
```
