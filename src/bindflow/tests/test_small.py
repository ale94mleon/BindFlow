#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Small-molecule parameterization through the external BindFlow-TOFF environment."""
import sys
import tempfile
from pathlib import Path

import pytest

from bindflow.utils import toff_env

tmp_dir = tempfile.TemporaryDirectory()
smi_file = Path(tmp_dir.name) / "ethane.smi"
smi_file.write_text("CC\n")

xfail_on_macos = pytest.mark.xfail(sys.platform == "darwin", reason="Expected to fail on macOS")


@pytest.fixture(scope="module")
def toff_prefix():
    """Locate the BindFlow-TOFF environment and check its TOFF version, or skip if it is missing."""
    try:
        prefix = toff_env.find_toff_prefix()
    except RuntimeError as e:
        pytest.skip(str(e))
    toff_env.check_toff_installation(prefix)
    return prefix


@pytest.mark.parametrize("ff_type, ff_code, resname", [
    ("openff", None, "OPE"),
    pytest.param("gaff", None, "GAF", marks=xfail_on_macos),
    pytest.param("espaloma", "espaloma-0.3.1", "ESP", marks=xfail_on_macos),
])
def test_parameterize(toff_prefix, ff_type, ff_code, resname):
    """Parameterize ethane with each force field type and check that the topology and coordinates are written."""
    if ff_type == "espaloma":
        try:
            toff_env.run_in_env(toff_prefix, ["python", "-c", "import espaloma"])
        except RuntimeError:
            pytest.skip("Espaloma is not installed in the TOFF environment")
    toff_env.parameterize(prefix=toff_prefix, input_mol=smi_file, out_dir=tmp_dir.name, mol_resi_name=resname,
                          force_field_type=ff_type, force_field_code=ff_code, hmr_factor=2.5)
    assert (Path(tmp_dir.name) / f"{resname}.top").is_file()
    assert (Path(tmp_dir.name) / f"{resname}.gro").is_file()
