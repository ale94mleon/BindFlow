"""FEP analysis: alchemlyb estimates, per-window sampling and the thermodynamic cycle.

The estimates use the GROMACS benzene data set of alchemtest (5 Coulomb and 16 vdW windows,
300 K). The reference values are regression values: they pin down the current results, so a
change in the subsampling or the estimators that moves them gets noticed.
"""
import json
import math
import warnings

import pytest

from bindflow.free_energy import fep_analysis

alchemtest_gmx = pytest.importorskip("alchemtest.gmx", reason="alchemtest provides the GROMACS benzene test data")

TEMPERATURE = 300  # K, temperature of the benzene simulations
# kcal/mol: (value, analytical error) for min_samples = 500
REFERENCE = {
    "Coulomb": {"MBAR": (1.8120, 0.0159), "TI": (1.8395, 0.0167)},
    "VDW": {"MBAR": (-1.7864, 0.0339), "TI": (-1.8262, 0.0372)},
}
SAMPLING_KEYS = {"window", "frames", "frame_interval_ps", "statistical_inefficiency", "statistical_inefficiency_ps",
                 "independent_samples", "step", "samples", "decorrelated"}


@pytest.fixture(scope="module")
def benzene():
    """Sorted xvg files of the benzene data set, per lambda type."""
    return {kind: sorted(files) for kind, files in alchemtest_gmx.load_benzene().data.items()}


@pytest.fixture(scope="module")
def benzene_results(benzene):
    """run_alchemlyb on each lambda type, computed once for the module."""
    return {kind: fep_analysis.run_alchemlyb(files, min_samples=500, temperature=TEMPERATURE)
            for kind, files in benzene.items()}


@pytest.mark.parametrize("kind", ["Coulomb", "VDW"])
@pytest.mark.parametrize("estimator", ["MBAR", "TI"])
def test_run_alchemlyb_benzene(benzene_results, kind, estimator):
    value, error = REFERENCE[kind][estimator]
    result = benzene_results[kind][estimator]
    assert result["value"] == pytest.approx(value, abs=1e-3)
    assert result["error"] == pytest.approx(error, abs=1e-3)


@pytest.mark.parametrize("kind", ["Coulomb", "VDW"])
def test_mbar_and_ti_agree(benzene_results, kind):
    assert abs(benzene_results[kind]["MBAR"]["value"] - benzene_results[kind]["TI"]["value"]) < 0.1


@pytest.mark.parametrize("kind", ["Coulomb", "VDW"])
def test_sampling_records(benzene, benzene_results, kind):
    sampling = benzene_results[kind]["sampling"]
    assert len(sampling) == len(benzene[kind])
    for window in sampling:
        assert set(window) == SAMPLING_KEYS
        assert window["samples"] == math.ceil(window["frames"] / window["step"])
        assert window["statistical_inefficiency_ps"] == pytest.approx(
            window["statistical_inefficiency"] * window["frame_interval_ps"])
        # The benzene windows are short-correlated: 500 samples never forces correlated frames
        assert window["decorrelated"]
    json.dumps(sampling)  # it is written to the contribution JSON files


def test_min_samples_floor_keeps_correlated_frames_and_warns(benzene):
    # Every benzene window has 4001 frames, so min_samples = 3000 forces step 1, below the
    # statistical inefficiency of most windows
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = fep_analysis.run_alchemlyb(benzene["Coulomb"], min_samples=3000, temperature=TEMPERATURE)
    assert all(window["step"] == 1 for window in result["sampling"])
    correlated = [window for window in result["sampling"] if not window["decorrelated"]]
    assert correlated
    assert any("kept correlated samples" in str(warning.message) for warning in caught)


def test_get_dG_contributions_writes_json(benzene, tmp_path):
    boresch = tmp_path / "dG_off.dat"
    boresch.write_text("7.5")
    out_json = tmp_path / "contributions.json"
    fep_analysis.get_dG_contributions(boresch_data=boresch, out_json_path=out_json, temperature=TEMPERATURE,
                                      vdw=benzene["VDW"], coul=benzene["Coulomb"])
    data = json.loads(out_json.read_text())
    assert set(data) == {"vdw", "coul", "boresch"}
    assert data["boresch"] == 7.5
    for kind in ("vdw", "coul"):
        assert set(data[kind]) == {"MBAR", "TI", "sampling"}


def test_get_dG_contributions_rejects_unknown_lambda_type(tmp_path):
    with pytest.raises(ValueError, match="not valid"):
        fep_analysis.get_dG_contributions(out_json_path=tmp_path / "x.json", restraints=["a.xvg"])


def _contribution(value, error):
    return {"MBAR": {"value": value, "error": error}, "TI": {"value": value + 0.1, "error": error}}


def test_get_dg_cycle(tmp_path):
    """ΔG = ligand(coul + vdw) − boresch + complex(vdw + coul) − complex(bonded); errors add in quadrature."""
    ligand = {"coul": _contribution(1.0, 0.1), "vdw": _contribution(2.0, 0.2)}
    complex_ = {"coul": _contribution(-10.0, 0.3), "vdw": _contribution(-20.0, 0.4),
                "bonded": _contribution(3.0, 0.5), "boresch": 7.0}
    (tmp_path / "ligand.json").write_text(json.dumps(ligand))
    (tmp_path / "complex.json").write_text(json.dumps(complex_))
    out_csv = tmp_path / "dG_results.csv"
    fep_analysis.get_dg_cycle(tmp_path / "ligand.json", tmp_path / "complex.json", out_csv)

    import pandas as pd
    result = pd.read_csv(out_csv, index_col=0)
    expected_error = math.sqrt(0.1**2 + 0.2**2 + 0.3**2 + 0.4**2 + 0.5**2)
    assert result.loc["value", "MBAR"] == pytest.approx(1.0 + 2.0 - 7.0 - 20.0 - 10.0 - 3.0)
    # TI: every contribution is 0.1 larger; bonded enters with a minus sign, so 4 × 0.1 − 0.1
    assert result.loc["value", "TI"] == pytest.approx(1.0 + 2.0 - 7.0 - 20.0 - 10.0 - 3.0 + 0.3)
    assert result.loc["std_dev", "MBAR"] == pytest.approx(expected_error)
    assert result.loc["std_dev", "TI"] == pytest.approx(expected_error)
