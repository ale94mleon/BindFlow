#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import tempfile
import sys
import pytest
from rdkit import Chem
from toff import Parameterize

tmp_dir = tempfile.TemporaryDirectory()
mol = Chem.MolFromSmiles('CC')


@pytest.mark.filterwarnings("ignore")
def test_Parameterize_openff():
    parameterizer = Parameterize(overwrite=True, out_dir=tmp_dir.name, force_field_type='openff', hmr_factor=2.5)
    parameterizer(input_mol=mol, mol_resi_name='OPE',)


@pytest.mark.xfail(sys.platform == "darwin", reason="This test is expected to fail on macOS")
@pytest.mark.filterwarnings("ignore")
def test_Parameterize_gaff():
    parameterizer = Parameterize(overwrite=True, out_dir=tmp_dir.name, force_field_type='gaff')
    parameterizer(input_mol=mol, mol_resi_name='GAF')

@pytest.mark.xfail(sys.platform == "darwin", reason="This test is expected to fail on macOS")
@pytest.mark.filterwarnings("ignore")
def test_Parameterize_espaloma():
    parameterizer = Parameterize(overwrite=True, out_dir=tmp_dir.name, force_field_type='espaloma', force_field_code='espaloma-0.3.1')
    parameterizer(input_mol=mol, mol_resi_name='ESP')


if __name__ == '__main__':
    pass
