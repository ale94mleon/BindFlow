import GMXMMPBSA.API
import math
import pandas as pd


def prettify_df(full_df):
    prettified_dict = {
        "name": [],

        "dG_pb_c2": [],
        "dG_pb_c2_err": [],
        "dG_pb_ie": [],
        "dG_pb_ie_err": [],
        "dG_pb_qh": [],
        "dG_pb_qh_err": [],

        "dG_gb_c2": [],
        "dG_gb_c2_err": [],
        "dG_gb_ie": [],
        "dG_gb_ie_err": [],
        "dG_gb_qh": [],
        "dG_gb_qh_err": [],
        }
    groups = full_df.groupby("name")
    for group_name, group_df in groups:
        prettified_dict["name"].append(group_name)
        prettified_dict["dG_pb_c2"].append(group_df["pb_c2_val"].mean())
        prettified_dict["dG_pb_c2_err"].append(group_df["pb_c2_val"].std())
        prettified_dict["dG_pb_ie"].append(group_df["pb_ie_val"].mean())
        prettified_dict["dG_pb_ie_err"].append(group_df["pb_ie_val"].std())
        prettified_dict["dG_pb_qh"].append(group_df["pb_qh_val"].mean())
        prettified_dict["dG_pb_qh_err"].append(group_df["pb_qh_val"].std())
        prettified_dict["dG_gb_c2"].append(group_df["gb_c2_val"].mean())
        prettified_dict["dG_gb_c2_err"].append(group_df["gb_c2_val"].std())
        prettified_dict["dG_gb_ie"].append(group_df["gb_ie_val"].mean())
        prettified_dict["dG_gb_ie_err"].append(group_df["gb_ie_val"].std())
        prettified_dict["dG_gb_qh"].append(group_df["gb_qh_val"].mean())
        prettified_dict["dG_gb_qh_err"].append(group_df["gb_qh_val"].std())
    
    return pd.DataFrame(prettified_dict)


def convert_format_flatten(df, ligand_name, replica, sample):
    res = { "name":[ligand_name], "replica":[replica], "sample":[sample], 
            "pb_c2_val":[], "pb_c2_err":[], "pb_ie_val":[], "pb_ie_err":[], "pb_qh_val":[], "pb_qh_err":[], 
            "gb_c2_val":[], "gb_c2_err":[], "gb_ie_val":[], "gb_ie_err":[], "gb_qh_val":[], "gb_qh_err":[]}
    
    def __set_row(name):
        extracted_row = df[df["method"] == name]
        res[name+"_val"].append(extracted_row.iat[0,1])
        res[name+"_err"].append(extracted_row.iat[0,2])
    __set_row("pb_c2")
    __set_row("pb_ie")
    __set_row("pb_qh")
    __set_row("gb_c2")
    __set_row("gb_ie")
    __set_row("gb_qh")
    
    return pd.DataFrame.from_dict(res)


def compute_dg(energy_value, error_energy, entropy_value, error_entropy):
    if energy_value is None or entropy_value is None:
        return None, None
    delta_g = energy_value+entropy_value
    if not error_energy is None and not error_entropy is None:
        delta_g_error = math.sqrt(error_energy**2 + error_entropy**2)
    else:
        delta_g_error = None
    return delta_g, delta_g_error


class GmxMmxbsaDataRetriever:
    def __init__(self, binary_api_file):
        """This class extracts the data generated from the GMX_MMPBSA program.
        Note that this class requires the binary data file generated after setting 
        the keep_files=0 or keep_files=1 parameter (keep_files=2 does not generate
        the required file)"""
        self.gmx_api = GMXMMPBSA.API.MMPBSA_API()
        self.gmx_api.setting_time()
        self.gmx_api.load_file(binary_api_file)

        self._entropies = self.__extract_entropies()
        self._energies = self.__extract_energies()
        self._delta_gs = self.__initialize_delta_g()

    
    def __extract_entropies(self):
        entropies = {"method": [], "value": [], "error": []}
        def __append_empty_entropy(name):
                entropies["method"].append(name)
                entropies["value"].append(None)
                entropies["error"].append(None)

        if "c2" in self.gmx_api.data["normal"].keys():
            if "pb" in self.gmx_api.data["normal"]["c2"].keys():
                entropies["method"].append("c2_pb")
                entropies["value"].append(self.gmx_api.data["normal"]["c2"]["pb"]["c2data"])
                entropies["error"].append(self.gmx_api.data["normal"]["c2"]["pb"]["c2_std"])
            else:
                __append_empty_entropy("c2_pb")
            if "gb" in self.gmx_api.data["normal"]["c2"].keys():
                entropies["method"].append("c2_gb")
                entropies["value"].append(self.gmx_api.data["normal"]["c2"]["gb"]["c2data"])
                entropies["error"].append(self.gmx_api.data["normal"]["c2"]["gb"]["c2_std"])
            else:
                __append_empty_entropy("c2_gb")
        else:
            __append_empty_entropy("c2_pb")
            __append_empty_entropy("c2_gb")
        
        if "ie" in self.gmx_api.data["normal"].keys():
            if "pb" in self.gmx_api.data["normal"]["ie"].keys():
                entropies["method"].append("ie_pb")
                entropies["value"].append(float(self.gmx_api.data["normal"]["ie"]["pb"]["iedata"].mean()))
                entropies["error"].append(float(self.gmx_api.data["normal"]["ie"]["pb"]["iedata"].std()))
            else:
                __append_empty_entropy("ie_pb")
            if "gb" in self.gmx_api.data["normal"]["ie"].keys():
                entropies["method"].append("ie_gb")
                entropies["value"].append(float(self.gmx_api.data["normal"]["ie"]["gb"]["iedata"].mean()))
                entropies["error"].append(float(self.gmx_api.data["normal"]["ie"]["gb"]["iedata"].std()))
            else:
                __append_empty_entropy("ie_gb")
        else:
            __append_empty_entropy("ie_pb")
            __append_empty_entropy("ie_gb")
        if "qh" in self.gmx_api.data["normal"].keys():
            entropies["method"].append("qh")
            entropies["value"].append(self.gmx_api.data["normal"]["qh"]["delta"]["TOTAL"])
            entropies["error"].append(None)
        else:
            __append_empty_entropy("qh")
        return entropies


    def __extract_energies(self):
        energies = {"method": [], "value": [], "error": []}
        def __append_empty_energies(name):
            energies["method"].append(name)
            energies["value"].append(None)
            energies["error"].append(None)
        
        if "pb" in self.gmx_api.data["normal"].keys():
            energies["method"].append("pb")
            energies["value"].append(float(self.gmx_api.data["normal"]["pb"]["delta"]["TOTAL"].mean()))
            energies["error"].append(float(self.gmx_api.data["normal"]["pb"]["delta"]["TOTAL"].std()))
        else:
            __append_empty_energies("pb")
        if "gb" in self.gmx_api.data["normal"].keys():
            energies["method"].append("gb")
            energies["value"].append(float(self.gmx_api.data["normal"]["gb"]["delta"]["TOTAL"].mean()))
            energies["error"].append(float(self.gmx_api.data["normal"]["gb"]["delta"]["TOTAL"].std()))
        else:
            __append_empty_energies("gb")
        return energies


    def __initialize_delta_g(self):
        delta_gs = {"method": [], "delta_g": [], "delta_g_error": []}
        df_entropy = pd.DataFrame.from_dict(self._entropies)
        df_energy = pd.DataFrame.from_dict(self._energies)
        
        def get_elements_from_energy(method):
            return [df_energy[df_energy["method"] == method]["value"].iat[0], df_energy[df_energy["method"] == method]["error"].iat[0]]
        def get_elements_from_entropy(method):
            return [df_entropy[df_entropy["method"] == method]["value"].iat[0], df_entropy[df_entropy["method"] == method]["error"].iat[0]]
        
        def create_dg_dict(energy_name, entropy_name, out_name):
            delta_g_val, delta_g_err = compute_dg(*get_elements_from_energy(energy_name),*get_elements_from_entropy(entropy_name))
            delta_gs["method"].append(out_name)
            delta_gs["delta_g"].append(delta_g_val)
            delta_gs["delta_g_error"].append(delta_g_err)

        create_dg_dict("pb", "c2_pb", "pb_c2")
        create_dg_dict("pb", "ie_pb", "pb_ie")
        create_dg_dict("pb", "qh", "pb_qh")
        create_dg_dict("gb", "c2_gb", "gb_c2")
        create_dg_dict("gb", "ie_gb", "gb_ie")
        create_dg_dict("gb", "qh", "gb_qh")

        return delta_gs
    
    def get_dg(self):
        return pd.DataFrame.from_dict(self._delta_gs)