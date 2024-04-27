import bindflow.utils.tools
import os
import json
import shutil
import GMXMMPBSA.API
import pandas as pd
import math

approach_path = config["out_approach_path"]

rule run_gmxmmpbsa:
    input:
        finished = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.finished",
        xtc = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod_noPBC.xtc",
        top = approach_path + "/{ligand_name}/input/complex/complex.top",
        mmpbsa_in = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/mmpbsa.in",
        ndx = approach_path + "/{ligand_name}/input/complex/index.ndx",
    output:
        gmxmmpbsa_res = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/COMPACT_MMXSA_RESULTS.mmxsa",
    params:
        in_tpr = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/prod.tpr",
    run:
        out_ligand_path = os.path.join(out_approach_path, wildcards.ligand_name)
        builder_dir = os.path.join(out_ligand_path, "builder")
        
        # TODO: Use a temporal directory instead of the builder_dir and pass all the files to where they are needed.
        # TODO: try to avoid changing directory, check how to handled output of gmx_MMPBSA
        cwd = os.getcwd()
        os.makedirs(builder_dir, exist_ok=True)
        os.chdir(builder_dir)
        try:
            # The index file generated in bindflow.preparation.system_builder.MakeInputs.__call__
            # will always have as first group receptor and as second group ligand
            # therefore, we can pass to the flag -cg <Receptor group> <Ligand group>" = -cg 1 2
            gmx_mmpbsa_command = f"gmx_MMPBSA -O -i {input.mmpbsa_in} -cs {params.in_tpr} -ci {input.ndx} -cg 1 2 -ct {input.xtc} -cp {input.top} -o res.dat -nogui"
            bindflow.utils.tools.run(gmx_mmpbsa_command)
            print(os.path.join(builder_dir,"COMPACT_MMXSA_RESULTS.mmxsa"))
            print(output.gmxmmpbsa_res)
            shutil.move(os.path.join(builder_dir,"COMPACT_MMXSA_RESULTS.mmxsa"), output.gmxmmpbsa_res)
        finally:
            os.chdir(cwd)
            try:
                shutil.rmtree(builder_dir)
            except FileNotFoundError:
                pass


rule mmpbsa_get_dg_cycle:
    input:
        gmxmmpbsa_res = approach_path + "/{ligand_name}/{replica}/complex/equil-mdsim/prod/COMPACT_MMXSA_RESULTS.mmxsa",
    output:
        out_file_path = approach_path + "/{ligand_name}/{replica}/dG_results.csv",
    run:
        gmx_api = GMXMMPBSA.API.MMPBSA_API()
        gmx_api.setting_time()
        gmx_api.load_file(input.gmxmmpbsa_res)
        
        entropies = {"method": [], "value": [], "error": []}
        def __apend_empty(name):
            entropies["method"].append(name)
            entropies["value"].append(None)
            entropies["error"].append(None)
        
        if "c2" in gmx_api.data["normal"].keys():
            if "pb" in gmx_api.data["normal"]["c2"].keys():
                entropies["method"].append("c2_pb")
                entropies["value"].append(gmx_api.data["normal"]["c2"]["pb"]["c2data"])
                entropies["error"].append(gmx_api.data["normal"]["c2"]["pb"]["c2_std"])
            else:
                __apend_empty("c2_pb")
            if "gb" in gmx_api.data["normal"]["c2"].keys():
                entropies["method"].append("c2_gb")
                entropies["value"].append(gmx_api.data["normal"]["c2"]["gb"]["c2data"])
                entropies["error"].append(gmx_api.data["normal"]["c2"]["gb"]["c2_std"])
            else:
                __apend_empty("c2_gb")
        else:
            __apend_empty("c2_pb")
            __apend_empty("c2_gb")
        
        if "ie" in gmx_api.data["normal"].keys():
            if "pb" in gmx_api.data["normal"]["ie"].keys():
                entropies["method"].append("ie_pb")
                entropies["value"].append(float(gmx_api.data["normal"]["ie"]["pb"]["iedata"].mean()))
                entropies["error"].append(float(gmx_api.data["normal"]["ie"]["pb"]["iedata"].std()))
            else:
                __apend_empty("ie_pb")
            if "gb" in gmx_api.data["normal"]["ie"].keys():
                entropies["method"].append("ie_gb")
                entropies["value"].append(float(gmx_api.data["normal"]["ie"]["gb"]["iedata"].mean()))
                entropies["error"].append(float(gmx_api.data["normal"]["ie"]["gb"]["iedata"].std()))
            else:
                __apend_empty("ie_gb")
        else:
            __apend_empty("ie_pb")
            __apend_empty("ie_gb")
        if "qh" in gmx_api.data["normal"].keys():
            entropies["method"].append("qh")
            entropies["value"].append(gmx_api.data["normal"]["qh"]["delta"]["TOTAL"])
            entropies["error"].append(None)
        else:
            __apend_empty("qh")

        energies = {"method": [], "value": [], "error": []}
        def __apend_empty_en(name):
            energies["method"].append(name)
            energies["value"].append(None)
            energies["error"].append(None)
        if "pb" in gmx_api.data["normal"].keys():
            energies["method"].append("pb")
            energies["value"].append(float(gmx_api.data["normal"]["pb"]["delta"]["TOTAL"].mean()))
            energies["error"].append(float(gmx_api.data["normal"]["pb"]["delta"]["TOTAL"].std()))
        else:
            __apend_empty_en("pb")
        if "gb" in gmx_api.data["normal"].keys():
            energies["method"].append("gb")
            energies["value"].append(float(gmx_api.data["normal"]["gb"]["delta"]["TOTAL"].mean()))
            energies["error"].append(float(gmx_api.data["normal"]["gb"]["delta"]["TOTAL"].std()))
        else:
            __apend_empty_en("gb")

        def compute_dg(energy_value, error_energy, entropy_value, error_entropy):
            if energy_value is None or entropy_value is None:
                return None, None
            delta_g = energy_value+entropy_value
            if not error_energy is None and not error_entropy is None:
                delta_g_error = math.sqrt(error_energy**2 + error_entropy**2)
            else:
                delta_g_error = None
            return delta_g, delta_g_error
        df_entropy = pd.DataFrame.from_dict(entropies)
        df_energy = pd.DataFrame.from_dict(energies)
        def get_elements_from_energy(method):
            return [df_energy[df_energy["method"] == method]["value"].iat[0], df_energy[df_energy["method"] == method]["error"].iat[0]]
        def get_elements_from_entropy(method):
            return [df_entropy[df_entropy["method"] == method]["value"].iat[0], df_entropy[df_entropy["method"] == method]["error"].iat[0]]
        delta_gs = {"method": [], "delta_g": [], "delta_g_error": []}
        def __create_dg_dict(energy_name, entropy_name, out_name):
            delta_g_val, delta_g_err = compute_dg(*get_elements_from_energy(energy_name),*get_elements_from_entropy(entropy_name))
            delta_gs["method"].append(out_name)
            delta_gs["delta_g"].append(delta_g_val)
            delta_gs["delta_g_error"].append(delta_g_err)
        __create_dg_dict("pb", "c2_pb", "pb_c2")
        __create_dg_dict("pb", "ie_pb", "pb_ie")
        __create_dg_dict("pb", "qh", "pb_qh")
        __create_dg_dict("gb", "c2_gb", "gb_c2")
        __create_dg_dict("gb", "ie_gb", "gb_ie")
        __create_dg_dict("gb", "qh", "gb_qh")
        pd.DataFrame.from_dict(delta_gs).to_csv(output.out_file_path, index=False)