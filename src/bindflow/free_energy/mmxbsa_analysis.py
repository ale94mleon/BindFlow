import GMXMMPBSA.API
import math
import pandas as pd


def prettify_df(full_df):
    groups = full_df.groupby(["name", "replica"])
    pretty_dict = {"name":[], "replica": [], 
                   "dg_c2_pb":[],"dg_c2_pb_sem":[],"dg_c2_gb":[],"dg_c2_gb_sem":[],"dg_ie_pb":[],"dg_ie_pb_sem":[],
                   "dg_ie_gb":[],"dg_ie_gb_sem":[],"dg_en_pb":[],"dg_en_pb_sem":[],"dg_en_gb":[],"dg_en_gb_sem":[],
                   "c2_pb":[],"c2_gb":[],"ie_pb":[],"ie_gb":[]}

    for group_keys, group_df in groups:
        group_name, group_replica = group_keys
        pretty_dict["name"].append(group_name)
        pretty_dict["replica"].append(group_replica)

        pretty_dict["dg_c2_pb"].append(group_df["dg_c2_pb"].mean())
        pretty_dict["dg_c2_pb_sem"].append(group_df["dg_c2_pb"].sem())
        pretty_dict["dg_c2_gb"].append(group_df["dg_c2_gb"].mean())
        pretty_dict["dg_c2_gb_sem"].append(group_df["dg_c2_gb"].sem())
        pretty_dict["dg_ie_pb"].append(group_df["dg_ie_pb"].mean())
        pretty_dict["dg_ie_pb_sem"].append(group_df["dg_ie_pb"].sem())
        pretty_dict["dg_ie_gb"].append(group_df["dg_ie_gb"].mean())
        pretty_dict["dg_ie_gb_sem"].append(group_df["dg_ie_gb"].sem())
        pretty_dict["dg_en_pb"].append(group_df["dg_en_pb"].mean())
        pretty_dict["dg_en_pb_sem"].append(group_df["dg_en_pb"].sem())
        pretty_dict["dg_en_gb"].append(group_df["dg_en_gb"].mean())
        pretty_dict["dg_en_gb_sem"].append(group_df["dg_en_gb"].sem())

        pretty_dict["c2_pb"].append(group_df["c2_pb"].mean())
        pretty_dict["c2_gb"].append(group_df["c2_gb"].mean())
        pretty_dict["ie_pb"].append(group_df["ie_pb"].mean())
        pretty_dict["ie_gb"].append(group_df["ie_gb"].mean())
    
    return pd.DataFrame(pretty_dict)


def convert_format_flatten(df, ligand_name, replica, sample):
    res = {"name":[ligand_name], "replica":[replica], "sample":[sample], 
           "dg_c2_pb": df["dg_c2_pb"], "dg_c2_gb": df["dg_c2_gb"], "dg_ie_pb": df["dg_ie_pb"],
           "dg_ie_gb": df["dg_ie_gb"], "dg_en_pb": df["dg_en_pb"], "dg_en_gb": df["dg_en_gb"], 
           "c2_pb": df["c2_pb"], "c2_gb": df["c2_gb"], "ie_pb": df["ie_pb"], "ie_gb": df["ie_gb"]}  
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

        self.__extract_entropies()
        self.__extract_energies()

    
    def __extract_entropies(self):
        if "c2" in self.gmx_api.data["normal"].keys():
            if "pb" in self.gmx_api.data["normal"]["c2"].keys():
                self.c2_pb = self.gmx_api.data["normal"]["c2"]["pb"]["c2data"]
            else:
                self.c2_pb = None
            if "gb" in self.gmx_api.data["normal"]["c2"].keys():
                self.c2_gb = self.gmx_api.data["normal"]["c2"]["gb"]["c2data"]
            else:
                self.c2_gb = None
        else:
            self.c2_pb = None
            self.c2_gb = None


        if "ie" in self.gmx_api.data["normal"].keys():
            if "pb" in self.gmx_api.data["normal"]["ie"].keys():
                self.ie_pb = float(self.gmx_api.data["normal"]["ie"]["pb"]["iedata"].mean())
            else:
                self.ie_pb = None
            if "gb" in self.gmx_api.data["normal"]["ie"].keys():
                self.ie_gb = float(self.gmx_api.data["normal"]["ie"]["gb"]["iedata"].mean())
            else:
                self.ie_gb = None
        else:
            self.ie_pb = None
            self.ie_gb = None

        if "qh" in self.gmx_api.data["normal"].keys():
            self.qh = self.gmx_api.data["normal"]["qh"]["delta"]["TOTAL"]
        else:
            self.qh = None
        
        return self.c2_pb, self.c2_gb, self.ie_pb, self.ie_gb, self.qh


    def __extract_energies(self):
        if "pb" in self.gmx_api.data["normal"].keys():
            self.pb_en = self.gmx_api.data["normal"]["pb"]["delta"]["TOTAL"]
        else:
            self.pb_en = None
        if "gb" in self.gmx_api.data["normal"].keys():
            self.gb_en = self.gmx_api.data["normal"]["gb"]["delta"]["TOTAL"]
        else:
            self.gb_en = None
        
        return self.pb_en, self.gb_en
    

    def store_dg(self, output_file, run_dir):
        # storing pb energies of each frame 
        pd.DataFrame(self.pb_en, columns=["delta_g_pb"]).to_csv(f"{run_dir}pb_energy_frames", index=False)

        # storing gb energies of each frame 
        pd.DataFrame(self.gb_en, columns=["delta_g_gb"]).to_csv(f"{run_dir}gb_energy_frames", index=False)

        delta_g_dict = {
            "dg_c2_pb": [self.pb_en.mean()+self.c2_pb if not self.pb_en is None and not self.c2_pb is None else None], 
            "dg_c2_gb": [self.gb_en.mean()+self.c2_gb if not self.gb_en is None and not self.c2_gb is None else None], 
            "dg_ie_pb": [self.pb_en.mean()+self.ie_pb if not self.pb_en is None and not self.ie_pb is None else None], 
            "dg_ie_gb": [self.gb_en.mean()+self.ie_gb if not self.gb_en is None and not self.ie_gb is None else None],
            "dg_en_pb": [self.pb_en.mean()            if not self.pb_en is None else None], 
            "dg_en_gb": [self.gb_en.mean()            if not self.gb_en is None else None],
            "c2_pb": [self.c2_pb if not self.c2_pb is None else None], 
            "c2_gb": [self.c2_gb if not self.c2_gb is None else None],
            "ie_pb": [self.ie_pb if not self.ie_pb is None else None], 
            "ie_gb": [self.ie_gb if not self.ie_gb is None else None],
        }
        pd.DataFrame.from_dict(delta_g_dict).to_csv(output_file, index=False)