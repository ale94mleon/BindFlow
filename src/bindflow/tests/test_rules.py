"""Snakemake rules: dry runs that check which files a rule asks for, without running anything."""
import re
import shutil
import subprocess
from pathlib import Path

import pytest

from bindflow import rules

pytestmark = pytest.mark.skipif(shutil.which("snakemake") is None, reason="needs the snakemake command")

RULES_DIR = Path(rules.__file__).parent
# Different number of windows per lambda type, so a rule that mixes them up asks for files that do not exist
LAMBDAS = {"ligand": {"vdw": 3, "coul": 2}, "complex": {"vdw": 4, "coul": 2, "bonded": 5}}

SNAKEFILE = """
import os
from pathlib import Path
from bindflow.mdp import mdp
from bindflow.utils import tools
config = {{"lambdas": {lambdas}}}
out_approach_path = "out"
threads = 1
include: "{rules}/fep/ligand/ana.smk"
include: "{rules}/fep/complex/ana.smk"

rule all:
    input: "out/lig/1/{leg}/fep/ana/dg_{leg}_contributions.json"
"""


def missing_inputs(snakefile: Path) -> set:
    """Dry-run the workflow and return the input files Snakemake reports as missing."""
    command = [shutil.which("snakemake"), "-n", "--cores", "1", "--snakefile", str(snakefile)]
    process = subprocess.run(command, cwd=snakefile.parent, capture_output=True, text=True,  # nosec B603: fixed arguments
                             timeout=300)
    output = process.stdout + process.stderr
    assert "MissingInputException" in output, output[-2000:]
    return set(re.findall(r"^\s+(out/\S+)$", output, flags=re.MULTILINE))


@pytest.mark.parametrize("leg", ["ligand", "complex"])
def test_fep_analysis_asks_for_the_mdp_of_each_existing_window(tmp_path, leg):
    lambdas = {name: {kind: [0.0] * n for kind, n in kinds.items()} for name, kinds in LAMBDAS.items()}
    snakefile = tmp_path / "Snakefile"
    snakefile.write_text(SNAKEFILE.format(lambdas=lambdas, rules=RULES_DIR, leg=leg))

    mdps = {path for path in missing_inputs(snakefile) if path.endswith("prod.mdp")}
    expected = {f"out/lig/1/{leg}/fep/simulation/{kind}.{state}/prod/prod.mdp"
                for kind, n in LAMBDAS[leg].items() for state in range(n)}
    assert mdps == expected
