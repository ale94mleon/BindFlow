import bindflow.utils.tools
import os
import json
import shutil
import GMXMMPBSA.API
import pandas as pd
import math
import tempfile

approach_path = config["out_approach_path"]
samples = list(map(str, range(1,1 + config["samples"])))

rule get_mmxbsa_result:
    input:
        mmxbsa_dat = expand(approach_path + "/{ligand_name}/{replica}/complex/mmpbsa/simulation/rep.{sample}/mmxbsa.dat", sample = samples, allow_missing = True)
    output:
        dG_mmxbsa_results = out_approach_path + "/{ligand_name}/{replica}/dG_mmxbsa_results.csv"
    run:
        gmx_api = GMXMMPBSA.API.MMPBSA_API()
        gmx_api.setting_time()
        gmx_api.load_file(input.mmxbsa_dat)
        
        entropies = {"method": [], "value": [], "error": []}
        def __append_empty(name):
            entropies["method"].append(name)
            entropies["value"].append(None)
            entropies["error"].append(None)
        
        if "c2" in gmx_api.data["normal"].keys():
            if "pb" in gmx_api.data["normal"]["c2"].keys():
                entropies["method"].append("c2_pb")
                entropies["value"].append(gmx_api.data["normal"]["c2"]["pb"]["c2data"])
                entropies["error"].append(gmx_api.data["normal"]["c2"]["pb"]["c2_std"])
            else:
                __append_empty("c2_pb")
            if "gb" in gmx_api.data["normal"]["c2"].keys():
                entropies["method"].append("c2_gb")
                entropies["value"].append(gmx_api.data["normal"]["c2"]["gb"]["c2data"])
                entropies["error"].append(gmx_api.data["normal"]["c2"]["gb"]["c2_std"])
            else:
                __append_empty("c2_gb")
        else:
            __append_empty("c2_pb")
            __append_empty("c2_gb")
        
        if "ie" in gmx_api.data["normal"].keys():
            if "pb" in gmx_api.data["normal"]["ie"].keys():
                entropies["method"].append("ie_pb")
                entropies["value"].append(float(gmx_api.data["normal"]["ie"]["pb"]["iedata"].mean()))
                entropies["error"].append(float(gmx_api.data["normal"]["ie"]["pb"]["iedata"].std()))
            else:
                __append_empty("ie_pb")
            if "gb" in gmx_api.data["normal"]["ie"].keys():
                entropies["method"].append("ie_gb")
                entropies["value"].append(float(gmx_api.data["normal"]["ie"]["gb"]["iedata"].mean()))
                entropies["error"].append(float(gmx_api.data["normal"]["ie"]["gb"]["iedata"].std()))
            else:
                __append_empty("ie_gb")
        else:
            __append_empty("ie_pb")
            __append_empty("ie_gb")
        if "qh" in gmx_api.data["normal"].keys():
            entropies["method"].append("qh")
            entropies["value"].append(gmx_api.data["normal"]["qh"]["delta"]["TOTAL"])
            entropies["error"].append(None)
        else:
            __append_empty("qh")

        energies = {"method": [], "value": [], "error": []}
        def __append_empty_en(name):
            energies["method"].append(name)
            energies["value"].append(None)
            energies["error"].append(None)
        if "pb" in gmx_api.data["normal"].keys():
            energies["method"].append("pb")
            energies["value"].append(float(gmx_api.data["normal"]["pb"]["delta"]["TOTAL"].mean()))
            energies["error"].append(float(gmx_api.data["normal"]["pb"]["delta"]["TOTAL"].std()))
        else:
            __append_empty_en("pb")
        if "gb" in gmx_api.data["normal"].keys():
            energies["method"].append("gb")
            energies["value"].append(float(gmx_api.data["normal"]["gb"]["delta"]["TOTAL"].mean()))
            energies["error"].append(float(gmx_api.data["normal"]["gb"]["delta"]["TOTAL"].std()))
        else:
            __append_empty_en("gb")

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