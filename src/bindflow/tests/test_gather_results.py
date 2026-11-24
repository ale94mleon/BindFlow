"""Gathering FEP results over ligands and replicas, on a synthetic results tree."""
import json
import math

import pandas as pd
import pytest

from bindflow.free_energy import gather_results

# ligand -> replica -> (MBAR value, MBAR error); TI is set to MBAR + 0.5 with the same error
RESULTS = {
    "lig.1": {"1": (-8.0, 0.2), "2": (-8.4, 0.3), "3": (-6.5, 0.1)},
    "lig2": {"1": (-5.0, 0.4), "2": (-5.2, 0.2)},
}


@pytest.fixture()
def results_tree(tmp_path):
    """<root>/<ligand>/<replica>/dG_results.csv plus the contribution JSON files of each replica."""
    for ligand, replicas in RESULTS.items():
        for replica, (value, error) in replicas.items():
            folder = tmp_path / ligand / replica
            folder.mkdir(parents=True)
            pd.DataFrame({"MBAR": {"value": value, "std_dev": error},
                          "TI": {"value": value + 0.5, "std_dev": error}}).to_csv(folder / "dG_results.csv")
            contribution = {"MBAR": {"value": value, "error": error}, "TI": {"value": value, "error": error}}
            for leg, kinds in (("complex", ("vdw", "coul", "bonded")), ("ligand", ("vdw", "coul"))):
                ana = folder / leg / "fep" / "ana"
                ana.mkdir(parents=True)
                data = {kind: contribution for kind in kinds}
                if leg == "complex":
                    data["boresch"] = 7.0
                (ana / f"dg_{leg}_contributions.json").write_text(json.dumps(data))
    # A folder without results must be ignored
    (tmp_path / "lig3" / "1").mkdir(parents=True)
    return tmp_path


def test_get_all_fep_dgs(results_tree):
    df = gather_results.get_all_fep_dgs(results_tree).set_index("ligand")
    assert set(df.index) == {"lig.1", "lig2"}
    for ligand, replicas in RESULTS.items():
        values = [value for value, _ in replicas.values()]
        errors = [error for _, error in replicas.values()]
        n = len(values)
        assert df.loc[ligand, "MBAR"] == pytest.approx(sum(values) / n)
        assert df.loc[ligand, "TI"] == pytest.approx(sum(values) / n + 0.5)
        # Standard error of the mean between replicas
        assert df.loc[ligand, "MBAR_sem"] == pytest.approx(pd.Series(values).std(ddof=1) / math.sqrt(n))
        # Analytical errors propagated through the average
        assert df.loc[ligand, "MBAR_uncertainty_propagation"] == pytest.approx(
            math.sqrt(sum((error / n) ** 2 for error in errors)))
        assert df.loc[ligand, "MBAR_num_replicas"] == n


def test_get_all_fep_dgs_without_results(tmp_path):
    assert gather_results.get_all_fep_dgs(tmp_path).empty


def test_get_raw_fep_data_keeps_full_names(results_tree):
    df = gather_results.get_raw_fep_data(results_tree)
    assert len(df) == sum(len(replicas) for replicas in RESULTS.values())
    # Folder names are kept whole: "lig.1" must not become "lig"
    assert set(df["ligand"]) == {"lig.1", "lig2"}
    assert set(df["replica"]) == {"1", "2", "3"}
    assert (df["boresch"] == 7.0).all()
