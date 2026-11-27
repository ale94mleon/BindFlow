"""Workflow configuration written by approach_flow and runners.calculate."""
import json

from bindflow import runners
from bindflow.orchestration.flow_builder import approach_flow
from bindflow.orchestration.generate_scheduler import FrontEnd
from bindflow.utils import toff_env


def _minimal_global_config(tmp_path, toff_prefix=None):
    protein = tmp_path / "protein.pdb"
    ligand = tmp_path / "lig1.sdf"
    protein.write_text("END\n")
    ligand.write_text("dummy\n")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    config = {
        "calculation_type": "fep",
        "out_approach_path": str(out_dir),
        "inputs": {
            "protein": {"conf": str(protein), "top": None, "ff": {"code": "amber99sb-ildn"}},
            "ligands": [{"conf": str(ligand), "top": None, "ff": {"type": "openff", "code": None}}],
            "cofactor": None,
            "membrane": None,
        },
        "water_model": "amber/tip3p",
        "host_name": "Protein",
        "host_selection": "protein and name CA",
        "fix_protein": True,
        "solv_d": 1.5,
        "solv_bt": "dodecahedron",
        "solv_rmin": 1,
        "solv_ion_conc": 0.15,
        "cofactor_on_protein": True,
        "cofactor_selection": "resname COF",
        "ligand_names": ["lig1"],
        "replicas": 1,
        "hmr_factor": 2.5,
        "custom_ff_path": None,
        "threads": 2,
        "extra_directives": {"dependencies": [], "mdrun": {"ligand": {}, "complex": {}}},
        "dt_max": 0.004,
        "job_prefix": "",
        "scheduler_class": FrontEnd,
        "cluster": {"options": {"calculation": None}},
        "num_jobs": 2,
    }
    if toff_prefix is not None:
        config["toff_prefix"] = toff_prefix
    return config, out_dir


def test_approach_flow_writes_toff_prefix_to_snake_conf(tmp_path):
    custom_prefix = tmp_path / "custom-toff"
    config, out_dir = _minimal_global_config(tmp_path, toff_prefix=custom_prefix)

    approach_flow(config, submit=False)

    snake_conf = json.loads((out_dir / "snake_conf.json").read_text())
    assert snake_conf["toff_prefix"] == str(custom_prefix)


def test_approach_flow_defaults_toff_prefix_to_none_when_absent(tmp_path):
    config, out_dir = _minimal_global_config(tmp_path)

    approach_flow(config, submit=False)

    snake_conf = json.loads((out_dir / "snake_conf.json").read_text())
    assert snake_conf["toff_prefix"] is None


def test_calculate_propagates_toff_prefix_to_snake_conf(tmp_path, monkeypatch):
    protein = tmp_path / "protein.pdb"
    ligand = tmp_path / "lig1.sdf"
    protein.write_text("END\n")
    ligand.write_text("dummy\n")

    custom_prefix = tmp_path / "my-toff-env"
    (custom_prefix / "bin").mkdir(parents=True)
    (custom_prefix / "bin" / "parameterize").write_text("#!/bin/sh\n")

    monkeypatch.setattr(runners, "check_gromacs_installation", lambda: None)
    monkeypatch.setattr(runners, "check_workflow_dependencies", lambda *args, **kwargs: None)
    monkeypatch.setattr(toff_env, "check_toff_installation", lambda prefix: "0.3.0")

    out_dir = tmp_path / "run"
    runners.calculate(
        calculation_type="fep",
        protein=protein,
        ligands=[ligand],
        scheduler_class=FrontEnd,
        out_root_folder_path=out_dir,
        submit=False,
        toff_prefix=custom_prefix,
    )

    snake_conf = json.loads((out_dir / "snake_conf.json").read_text())
    assert snake_conf["toff_prefix"] == str(custom_prefix.resolve())
