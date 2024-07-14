from pathlib import Path

import GMXMMPBSA.API
import pandas as pd


def convert_format_flatten(df, ligand_name, replica, sample):
    res = {
        "name": [ligand_name],
        "replica": [int(replica)],
        "sample": [int(sample)],
        "dg_c2_pb": df["dg_c2_pb"],
        "dg_c2_gb": df["dg_c2_gb"],
        "dg_ie_pb": df["dg_ie_pb"],
        "dg_ie_gb": df["dg_ie_gb"],
        "dg_qh_pb": df["dg_qh_pb"],
        "dg_qh_gb": df["dg_qh_gb"],
        "dg_en_pb": df["dg_en_pb"],
        "dg_en_gb": df["dg_en_gb"],
        "c2_pb": df["c2_pb"],
        "c2_gb": df["c2_gb"],
        "ie_pb": df["ie_pb"],
        "ie_gb": df["ie_gb"],
        "qh": df["qh"],

        "gb_energy_complex_bond": df["gb_energy_complex_bond"],
        "gb_energy_complex_angle": df["gb_energy_complex_angle"],
        "gb_energy_complex_dihed": df["gb_energy_complex_dihed"],
        "gb_energy_complex_vdwaals": df["gb_energy_complex_vdwaals"],
        "gb_energy_complex_eel": df["gb_energy_complex_eel"],
        "gb_energy_complex_14vdw": df["gb_energy_complex_14vdw"],
        "gb_energy_complex_14eel": df["gb_energy_complex_14eel"],
        "gb_energy_complex_egb": df["gb_energy_complex_egb"],
        "gb_energy_complex_esurf": df["gb_energy_complex_esurf"],
        "gb_energy_complex_ggas": df["gb_energy_complex_ggas"],
        "gb_energy_complex_gsolv": df["gb_energy_complex_gsolv"],
        "gb_energy_complex_total": df["gb_energy_complex_total"],

        "gb_energy_receptor_bond": df["gb_energy_receptor_bond"],
        "gb_energy_receptor_angle": df["gb_energy_receptor_angle"],
        "gb_energy_receptor_dihed": df["gb_energy_receptor_dihed"],
        "gb_energy_receptor_vdwaals": df["gb_energy_receptor_vdwaals"],
        "gb_energy_receptor_eel": df["gb_energy_receptor_eel"],
        "gb_energy_receptor_14vdw": df["gb_energy_receptor_14vdw"],
        "gb_energy_receptor_14eel": df["gb_energy_receptor_14eel"],
        "gb_energy_receptor_egb": df["gb_energy_receptor_egb"],
        "gb_energy_receptor_esurf": df["gb_energy_receptor_esurf"],
        "gb_energy_receptor_ggas": df["gb_energy_receptor_ggas"],
        "gb_energy_receptor_gsolv": df["gb_energy_receptor_gsolv"],
        "gb_energy_receptor_total": df["gb_energy_receptor_total"],

        "gb_energy_ligand_bond": df["gb_energy_ligand_bond"],
        "gb_energy_ligand_angle": df["gb_energy_ligand_angle"],
        "gb_energy_ligand_dihed": df["gb_energy_ligand_dihed"],
        "gb_energy_ligand_vdwaals": df["gb_energy_ligand_vdwaals"],
        "gb_energy_ligand_eel": df["gb_energy_ligand_eel"],
        "gb_energy_ligand_14vdw": df["gb_energy_ligand_14vdw"],
        "gb_energy_ligand_14eel": df["gb_energy_ligand_14eel"],
        "gb_energy_ligand_egb": df["gb_energy_ligand_egb"],
        "gb_energy_ligand_esurf": df["gb_energy_ligand_esurf"],
        "gb_energy_ligand_ggas": df["gb_energy_ligand_ggas"],
        "gb_energy_ligand_gsolv": df["gb_energy_ligand_gsolv"],
        "gb_energy_ligand_total": df["gb_energy_ligand_total"],

        "gb_energy_delta_bond": df["gb_energy_delta_bond"],
        "gb_energy_delta_angle": df["gb_energy_delta_angle"],
        "gb_energy_delta_dihed": df["gb_energy_delta_dihed"],
        "gb_energy_delta_vdwaals": df["gb_energy_delta_vdwaals"],
        "gb_energy_delta_eel": df["gb_energy_delta_eel"],
        "gb_energy_delta_14vdw": df["gb_energy_delta_14vdw"],
        "gb_energy_delta_14eel": df["gb_energy_delta_14eel"],
        "gb_energy_delta_egb": df["gb_energy_delta_egb"],
        "gb_energy_delta_esurf": df["gb_energy_delta_esurf"],
        "gb_energy_delta_ggas": df["gb_energy_delta_ggas"],
        "gb_energy_delta_gsolv": df["gb_energy_delta_gsolv"],
        "gb_energy_delta_total": df["gb_energy_delta_total"],


        "pb_energy_complex_bond": df["pb_energy_complex_bond"],
        "pb_energy_complex_angle": df["pb_energy_complex_angle"],
        "pb_energy_complex_dihed": df["pb_energy_complex_dihed"],
        "pb_energy_complex_vdwaals": df["pb_energy_complex_vdwaals"],
        "pb_energy_complex_eel": df["pb_energy_complex_eel"],
        "pb_energy_complex_14vdw": df["pb_energy_complex_14vdw"],
        "pb_energy_complex_14eel": df["pb_energy_complex_14eel"],
        "pb_energy_complex_epb": df["pb_energy_complex_epb"],
        "pb_energy_complex_enpolar": df["pb_energy_complex_enpolar"],
        "pb_energy_complex_edisper": df["pb_energy_complex_edisper"],
        "pb_energy_complex_ggas": df["pb_energy_complex_ggas"],
        "pb_energy_complex_gsolv": df["pb_energy_complex_gsolv"],
        "pb_energy_complex_total": df["pb_energy_complex_total"],

        "pb_energy_receptor_bond": df["pb_energy_receptor_bond"],
        "pb_energy_receptor_angle": df["pb_energy_receptor_angle"],
        "pb_energy_receptor_dihed": df["pb_energy_receptor_dihed"],
        "pb_energy_receptor_vdwaals": df["pb_energy_receptor_vdwaals"],
        "pb_energy_receptor_eel": df["pb_energy_receptor_eel"],
        "pb_energy_receptor_14vdw": df["pb_energy_receptor_14vdw"],
        "pb_energy_receptor_14eel": df["pb_energy_receptor_14eel"],
        "pb_energy_receptor_epb": df["pb_energy_receptor_epb"],
        "pb_energy_receptor_enpolar": df["pb_energy_receptor_enpolar"],
        "pb_energy_receptor_edisper": df["pb_energy_receptor_edisper"],
        "pb_energy_receptor_ggas": df["pb_energy_receptor_ggas"],
        "pb_energy_receptor_gsolv": df["pb_energy_receptor_gsolv"],
        "pb_energy_receptor_total": df["pb_energy_receptor_total"],

        "pb_energy_ligand_bond": df["pb_energy_ligand_bond"],
        "pb_energy_ligand_angle": df["pb_energy_ligand_angle"],
        "pb_energy_ligand_dihed": df["pb_energy_ligand_dihed"],
        "pb_energy_ligand_vdwaals": df["pb_energy_ligand_vdwaals"],
        "pb_energy_ligand_eel": df["pb_energy_ligand_eel"],
        "pb_energy_ligand_14vdw": df["pb_energy_ligand_14vdw"],
        "pb_energy_ligand_14eel": df["pb_energy_ligand_14eel"],
        "pb_energy_ligand_epb": df["pb_energy_ligand_epb"],
        "pb_energy_ligand_enpolar": df["pb_energy_ligand_enpolar"],
        "pb_energy_ligand_edisper": df["pb_energy_ligand_edisper"],
        "pb_energy_ligand_ggas": df["pb_energy_ligand_ggas"],
        "pb_energy_ligand_gsolv": df["pb_energy_ligand_gsolv"],
        "pb_energy_ligand_total": df["pb_energy_ligand_total"],

        "pb_energy_delta_bond": df["pb_energy_delta_bond"],
        "pb_energy_delta_angle": df["pb_energy_delta_angle"],
        "pb_energy_delta_dihed": df["pb_energy_delta_dihed"],
        "pb_energy_delta_vdwaals": df["pb_energy_delta_vdwaals"],
        "pb_energy_delta_eel": df["pb_energy_delta_eel"],
        "pb_energy_delta_14vdw": df["pb_energy_delta_14vdw"],
        "pb_energy_delta_14eel": df["pb_energy_delta_14eel"],
        "pb_energy_delta_epb": df["pb_energy_delta_epb"],
        "pb_energy_delta_enpolar": df["pb_energy_delta_enpolar"],
        "pb_energy_delta_edisper": df["pb_energy_delta_edisper"],
        "pb_energy_delta_ggas": df["pb_energy_delta_ggas"],
        "pb_energy_delta_gsolv": df["pb_energy_delta_gsolv"],
        "pb_energy_delta_total": df["pb_energy_delta_total"],
    }
    return pd.DataFrame.from_dict(res)


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
        self.__extract_others()

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

    def __extract_others(self):
        if "gb" in self.gmx_api.data["normal"].keys():
            self.gb_energy_complex_bond = float(self.gmx_api.data["normal"]["gb"]["complex"]["BOND"].mean())
            self.gb_energy_complex_angle = float(self.gmx_api.data["normal"]["gb"]["complex"]["ANGLE"].mean())
            self.gb_energy_complex_dihed = float(self.gmx_api.data["normal"]["gb"]["complex"]["DIHED"].mean())
            self.gb_energy_complex_vdwaals = float(self.gmx_api.data["normal"]["gb"]["complex"]["VDWAALS"].mean())
            self.gb_energy_complex_eel = float(self.gmx_api.data["normal"]["gb"]["complex"]["EEL"].mean())
            self.gb_energy_complex_14vdw = float(self.gmx_api.data["normal"]["gb"]["complex"]["1-4 VDW"].mean())
            self.gb_energy_complex_14eel = float(self.gmx_api.data["normal"]["gb"]["complex"]["1-4 EEL"].mean())
            self.gb_energy_complex_egb = float(self.gmx_api.data["normal"]["gb"]["complex"]["EGB"].mean())
            self.gb_energy_complex_esurf = float(self.gmx_api.data["normal"]["gb"]["complex"]["ESURF"].mean())
            self.gb_energy_complex_ggas = float(self.gmx_api.data["normal"]["gb"]["complex"]["GGAS"].mean())
            self.gb_energy_complex_gsolv = float(self.gmx_api.data["normal"]["gb"]["complex"]["GSOLV"].mean())
            self.gb_energy_complex_total = float(self.gmx_api.data["normal"]["gb"]["complex"]["TOTAL"].mean())

            self.gb_energy_receptor_bond = float(self.gmx_api.data["normal"]["gb"]["receptor"]["BOND"].mean())
            self.gb_energy_receptor_angle = float(self.gmx_api.data["normal"]["gb"]["receptor"]["ANGLE"].mean())
            self.gb_energy_receptor_dihed = float(self.gmx_api.data["normal"]["gb"]["receptor"]["DIHED"].mean())
            self.gb_energy_receptor_vdwaals = float(self.gmx_api.data["normal"]["gb"]["receptor"]["VDWAALS"].mean())
            self.gb_energy_receptor_eel = float(self.gmx_api.data["normal"]["gb"]["receptor"]["EEL"].mean())
            self.gb_energy_receptor_14vdw = float(self.gmx_api.data["normal"]["gb"]["receptor"]["1-4 VDW"].mean())
            self.gb_energy_receptor_14eel = float(self.gmx_api.data["normal"]["gb"]["receptor"]["1-4 EEL"].mean())
            self.gb_energy_receptor_egb = float(self.gmx_api.data["normal"]["gb"]["receptor"]["EGB"].mean())
            self.gb_energy_receptor_esurf = float(self.gmx_api.data["normal"]["gb"]["receptor"]["ESURF"].mean())
            self.gb_energy_receptor_ggas = float(self.gmx_api.data["normal"]["gb"]["receptor"]["GGAS"].mean())
            self.gb_energy_receptor_gsolv = float(self.gmx_api.data["normal"]["gb"]["receptor"]["GSOLV"].mean())
            self.gb_energy_receptor_total = float(self.gmx_api.data["normal"]["gb"]["receptor"]["TOTAL"].mean())

            self.gb_energy_ligand_bond = float(self.gmx_api.data["normal"]["gb"]["ligand"]["BOND"].mean())
            self.gb_energy_ligand_angle = float(self.gmx_api.data["normal"]["gb"]["ligand"]["ANGLE"].mean())
            self.gb_energy_ligand_dihed = float(self.gmx_api.data["normal"]["gb"]["ligand"]["DIHED"].mean())
            self.gb_energy_ligand_vdwaals = float(self.gmx_api.data["normal"]["gb"]["ligand"]["VDWAALS"].mean())
            self.gb_energy_ligand_eel = float(self.gmx_api.data["normal"]["gb"]["ligand"]["EEL"].mean())
            self.gb_energy_ligand_14vdw = float(self.gmx_api.data["normal"]["gb"]["ligand"]["1-4 VDW"].mean())
            self.gb_energy_ligand_14eel = float(self.gmx_api.data["normal"]["gb"]["ligand"]["1-4 EEL"].mean())
            self.gb_energy_ligand_egb = float(self.gmx_api.data["normal"]["gb"]["ligand"]["EGB"].mean())
            self.gb_energy_ligand_esurf = float(self.gmx_api.data["normal"]["gb"]["ligand"]["ESURF"].mean())
            self.gb_energy_ligand_ggas = float(self.gmx_api.data["normal"]["gb"]["ligand"]["GGAS"].mean())
            self.gb_energy_ligand_gsolv = float(self.gmx_api.data["normal"]["gb"]["ligand"]["GSOLV"].mean())
            self.gb_energy_ligand_total = float(self.gmx_api.data["normal"]["gb"]["ligand"]["TOTAL"].mean())

            self.gb_energy_delta_bond = float(self.gmx_api.data["normal"]["gb"]["delta"]["BOND"].mean())
            self.gb_energy_delta_angle = float(self.gmx_api.data["normal"]["gb"]["delta"]["ANGLE"].mean())
            self.gb_energy_delta_dihed = float(self.gmx_api.data["normal"]["gb"]["delta"]["DIHED"].mean())
            self.gb_energy_delta_vdwaals = float(self.gmx_api.data["normal"]["gb"]["delta"]["VDWAALS"].mean())
            self.gb_energy_delta_eel = float(self.gmx_api.data["normal"]["gb"]["delta"]["EEL"].mean())
            self.gb_energy_delta_14vdw = float(self.gmx_api.data["normal"]["gb"]["delta"]["1-4 VDW"].mean())
            self.gb_energy_delta_14eel = float(self.gmx_api.data["normal"]["gb"]["delta"]["1-4 EEL"].mean())
            self.gb_energy_delta_egb = float(self.gmx_api.data["normal"]["gb"]["delta"]["EGB"].mean())
            self.gb_energy_delta_esurf = float(self.gmx_api.data["normal"]["gb"]["delta"]["ESURF"].mean())
            self.gb_energy_delta_ggas = float(self.gmx_api.data["normal"]["gb"]["delta"]["GGAS"].mean())
            self.gb_energy_delta_gsolv = float(self.gmx_api.data["normal"]["gb"]["delta"]["GSOLV"].mean())
            self.gb_energy_delta_total = float(self.gmx_api.data["normal"]["gb"]["delta"]["TOTAL"].mean())
        else:
            self.gb_energy_complex_bond = None
            self.gb_energy_complex_angle = None
            self.gb_energy_complex_dihed = None
            self.gb_energy_complex_vdwaals = None
            self.gb_energy_complex_eel = None
            self.gb_energy_complex_14vdw = None
            self.gb_energy_complex_14eel = None
            self.gb_energy_complex_egb = None
            self.gb_energy_complex_esurf = None
            self.gb_energy_complex_ggas = None
            self.gb_energy_complex_gsolv = None
            self.gb_energy_complex_total = None

            self.gb_energy_receptor_bond = None
            self.gb_energy_receptor_angle = None
            self.gb_energy_receptor_dihed = None
            self.gb_energy_receptor_vdwaals = None
            self.gb_energy_receptor_eel = None
            self.gb_energy_receptor_14vdw = None
            self.gb_energy_receptor_14eel = None
            self.gb_energy_receptor_egb = None
            self.gb_energy_receptor_esurf = None
            self.gb_energy_receptor_ggas = None
            self.gb_energy_receptor_gsolv = None
            self.gb_energy_receptor_total = None

            self.gb_energy_ligand_bond = None
            self.gb_energy_ligand_angle = None
            self.gb_energy_ligand_dihed = None
            self.gb_energy_ligand_vdwaals = None
            self.gb_energy_ligand_eel = None
            self.gb_energy_ligand_14vdw = None
            self.gb_energy_ligand_14eel = None
            self.gb_energy_ligand_egb = None
            self.gb_energy_ligand_esurf = None
            self.gb_energy_ligand_ggas = None
            self.gb_energy_ligand_gsolv = None
            self.gb_energy_ligand_total = None

            self.gb_energy_delta_bond = None
            self.gb_energy_delta_angle = None
            self.gb_energy_delta_dihed = None
            self.gb_energy_delta_vdwaals = None
            self.gb_energy_delta_eel = None
            self.gb_energy_delta_14vdw = None
            self.gb_energy_delta_14eel = None
            self.gb_energy_delta_egb = None
            self.gb_energy_delta_esurf = None
            self.gb_energy_delta_ggas = None
            self.gb_energy_delta_gsolv = None
            self.gb_energy_delta_total = None

        if "pb" in self.gmx_api.data["normal"].keys():
            self.pb_energy_complex_bond = float(self.gmx_api.data["normal"]["pb"]["complex"]["BOND"].mean())
            self.pb_energy_complex_angle = float(self.gmx_api.data["normal"]["pb"]["complex"]["ANGLE"].mean())
            self.pb_energy_complex_dihed = float(self.gmx_api.data["normal"]["pb"]["complex"]["DIHED"].mean())
            self.pb_energy_complex_vdwaals = float(self.gmx_api.data["normal"]["pb"]["complex"]["VDWAALS"].mean())
            self.pb_energy_complex_eel = float(self.gmx_api.data["normal"]["pb"]["complex"]["EEL"].mean())
            self.pb_energy_complex_14vdw = float(self.gmx_api.data["normal"]["pb"]["complex"]["1-4 VDW"].mean())
            self.pb_energy_complex_14eel = float(self.gmx_api.data["normal"]["pb"]["complex"]["1-4 EEL"].mean())
            self.pb_energy_complex_epb = float(self.gmx_api.data["normal"]["pb"]["complex"]["EPB"].mean())
            self.pb_energy_complex_enpolar = float(self.gmx_api.data["normal"]["pb"]["complex"]["ENPOLAR"].mean())
            self.pb_energy_complex_edisper = float(self.gmx_api.data["normal"]["pb"]["complex"]["EDISPER"].mean())
            self.pb_energy_complex_ggas = float(self.gmx_api.data["normal"]["pb"]["complex"]["GGAS"].mean())
            self.pb_energy_complex_gsolv = float(self.gmx_api.data["normal"]["pb"]["complex"]["GSOLV"].mean())
            self.pb_energy_complex_total = float(self.gmx_api.data["normal"]["pb"]["complex"]["TOTAL"].mean())

            self.pb_energy_receptor_bond = float(self.gmx_api.data["normal"]["pb"]["receptor"]["BOND"].mean())
            self.pb_energy_receptor_angle = float(self.gmx_api.data["normal"]["pb"]["receptor"]["ANGLE"].mean())
            self.pb_energy_receptor_dihed = float(self.gmx_api.data["normal"]["pb"]["receptor"]["DIHED"].mean())
            self.pb_energy_receptor_vdwaals = float(self.gmx_api.data["normal"]["pb"]["receptor"]["VDWAALS"].mean())
            self.pb_energy_receptor_eel = float(self.gmx_api.data["normal"]["pb"]["receptor"]["EEL"].mean())
            self.pb_energy_receptor_14vdw = float(self.gmx_api.data["normal"]["pb"]["receptor"]["1-4 VDW"].mean())
            self.pb_energy_receptor_14eel = float(self.gmx_api.data["normal"]["pb"]["receptor"]["1-4 EEL"].mean())
            self.pb_energy_receptor_epb = float(self.gmx_api.data["normal"]["pb"]["receptor"]["EPB"].mean())
            self.pb_energy_receptor_enpolar = float(self.gmx_api.data["normal"]["pb"]["receptor"]["ENPOLAR"].mean())
            self.pb_energy_receptor_edisper = float(self.gmx_api.data["normal"]["pb"]["receptor"]["EDISPER"].mean())
            self.pb_energy_receptor_ggas = float(self.gmx_api.data["normal"]["pb"]["receptor"]["GGAS"].mean())
            self.pb_energy_receptor_gsolv = float(self.gmx_api.data["normal"]["pb"]["receptor"]["GSOLV"].mean())
            self.pb_energy_receptor_total = float(self.gmx_api.data["normal"]["pb"]["receptor"]["TOTAL"].mean())

            self.pb_energy_ligand_bond = float(self.gmx_api.data["normal"]["pb"]["ligand"]["BOND"].mean())
            self.pb_energy_ligand_angle = float(self.gmx_api.data["normal"]["pb"]["ligand"]["ANGLE"].mean())
            self.pb_energy_ligand_dihed = float(self.gmx_api.data["normal"]["pb"]["ligand"]["DIHED"].mean())
            self.pb_energy_ligand_vdwaals = float(self.gmx_api.data["normal"]["pb"]["ligand"]["VDWAALS"].mean())
            self.pb_energy_ligand_eel = float(self.gmx_api.data["normal"]["pb"]["ligand"]["EEL"].mean())
            self.pb_energy_ligand_14vdw = float(self.gmx_api.data["normal"]["pb"]["ligand"]["1-4 VDW"].mean())
            self.pb_energy_ligand_14eel = float(self.gmx_api.data["normal"]["pb"]["ligand"]["1-4 EEL"].mean())
            self.pb_energy_ligand_epb = float(self.gmx_api.data["normal"]["pb"]["ligand"]["EPB"].mean())
            self.pb_energy_ligand_enpolar = float(self.gmx_api.data["normal"]["pb"]["ligand"]["ENPOLAR"].mean())
            self.pb_energy_ligand_edisper = float(self.gmx_api.data["normal"]["pb"]["ligand"]["EDISPER"].mean())
            self.pb_energy_ligand_ggas = float(self.gmx_api.data["normal"]["pb"]["ligand"]["GGAS"].mean())
            self.pb_energy_ligand_gsolv = float(self.gmx_api.data["normal"]["pb"]["ligand"]["GSOLV"].mean())
            self.pb_energy_ligand_total = float(self.gmx_api.data["normal"]["pb"]["ligand"]["TOTAL"].mean())

            self.pb_energy_delta_bond = float(self.gmx_api.data["normal"]["pb"]["delta"]["BOND"].mean())
            self.pb_energy_delta_angle = float(self.gmx_api.data["normal"]["pb"]["delta"]["ANGLE"].mean())
            self.pb_energy_delta_dihed = float(self.gmx_api.data["normal"]["pb"]["delta"]["DIHED"].mean())
            self.pb_energy_delta_vdwaals = float(self.gmx_api.data["normal"]["pb"]["delta"]["VDWAALS"].mean())
            self.pb_energy_delta_eel = float(self.gmx_api.data["normal"]["pb"]["delta"]["EEL"].mean())
            self.pb_energy_delta_14vdw = float(self.gmx_api.data["normal"]["pb"]["delta"]["1-4 VDW"].mean())
            self.pb_energy_delta_14eel = float(self.gmx_api.data["normal"]["pb"]["delta"]["1-4 EEL"].mean())
            self.pb_energy_delta_epb = float(self.gmx_api.data["normal"]["pb"]["delta"]["EPB"].mean())
            self.pb_energy_delta_enpolar = float(self.gmx_api.data["normal"]["pb"]["delta"]["ENPOLAR"].mean())
            self.pb_energy_delta_edisper = float(self.gmx_api.data["normal"]["pb"]["delta"]["EDISPER"].mean())
            self.pb_energy_delta_ggas = float(self.gmx_api.data["normal"]["pb"]["delta"]["GGAS"].mean())
            self.pb_energy_delta_gsolv = float(self.gmx_api.data["normal"]["pb"]["delta"]["GSOLV"].mean())
            self.pb_energy_delta_total = float(self.gmx_api.data["normal"]["pb"]["delta"]["TOTAL"].mean())
        else:
            self.pb_energy_complex_bond = None
            self.pb_energy_complex_angle = None
            self.pb_energy_complex_dihed = None
            self.pb_energy_complex_vdwaals = None
            self.pb_energy_complex_eel = None
            self.pb_energy_complex_14vdw = None
            self.pb_energy_complex_14eel = None
            self.pb_energy_complex_epb = None
            self.pb_energy_complex_enpolar = None
            self.pb_energy_complex_edisper = None
            self.pb_energy_complex_ggas = None
            self.pb_energy_complex_gsolv = None
            self.pb_energy_complex_total = None

            self.pb_energy_receptor_bond = None
            self.pb_energy_receptor_angle = None
            self.pb_energy_receptor_dihed = None
            self.pb_energy_receptor_vdwaals = None
            self.pb_energy_receptor_eel = None
            self.pb_energy_receptor_14vdw = None
            self.pb_energy_receptor_14eel = None
            self.pb_energy_receptor_epb = None
            self.pb_energy_receptor_enpolar = None
            self.pb_energy_receptor_edisper = None
            self.pb_energy_receptor_ggas = None
            self.pb_energy_receptor_gsolv = None
            self.pb_energy_receptor_total = None

            self.pb_energy_ligand_bond = None
            self.pb_energy_ligand_angle = None
            self.pb_energy_ligand_dihed = None
            self.pb_energy_ligand_vdwaals = None
            self.pb_energy_ligand_eel = None
            self.pb_energy_ligand_14vdw = None
            self.pb_energy_ligand_14eel = None
            self.pb_energy_ligand_epb = None
            self.pb_energy_ligand_enpolar = None
            self.pb_energy_ligand_edisper = None
            self.pb_energy_ligand_ggas = None
            self.pb_energy_ligand_gsolv = None
            self.pb_energy_ligand_total = None

            self.pb_energy_delta_bond = None
            self.pb_energy_delta_angle = None
            self.pb_energy_delta_dihed = None
            self.pb_energy_delta_vdwaals = None
            self.pb_energy_delta_eel = None
            self.pb_energy_delta_14vdw = None
            self.pb_energy_delta_14eel = None
            self.pb_energy_delta_epb = None
            self.pb_energy_delta_enpolar = None
            self.pb_energy_delta_edisper = None
            self.pb_energy_delta_ggas = None
            self.pb_energy_delta_gsolv = None
            self.pb_energy_delta_total = None

    def store_dg(self, output_file, run_dir):
        # storing pb energies of each frame
        pd.DataFrame(self.pb_en, columns=["delta_g_pb"]).to_csv(Path(run_dir)/"pb_energy_frames.csv", index=False)

        # storing gb energies of each frame
        pd.DataFrame(self.gb_en, columns=["delta_g_gb"]).to_csv(Path(run_dir)/"gb_energy_frames.csv", index=False)

        delta_g_dict = {
            "dg_c2_pb": [self.pb_en.mean() + self.c2_pb if (self.pb_en is not None and self.c2_pb is not None) else None],
            "dg_c2_gb": [self.gb_en.mean() + self.c2_gb if (self.gb_en is not None and self.c2_gb is not None) else None],
            "dg_ie_pb": [self.pb_en.mean() + self.ie_pb if (self.pb_en is not None and self.ie_pb is not None) else None],
            "dg_ie_gb": [self.gb_en.mean() + self.ie_gb if (self.gb_en is not None and self.ie_gb is not None) else None],
            "dg_qh_pb": [self.pb_en.mean() + self.qh if (self.pb_en is not None and self.qh is not None) else None],
            "dg_qh_gb": [self.gb_en.mean() + self.qh if (self.gb_en is not None and self.qh is not None) else None],
            "dg_en_pb": [self.pb_en.mean() if self.pb_en is not None else None],
            "dg_en_gb": [self.gb_en.mean() if self.gb_en is not None else None],
            "c2_pb": [self.c2_pb if self.c2_pb is not None else None],
            "c2_gb": [self.c2_gb if self.c2_gb is not None else None],
            "ie_pb": [self.ie_pb if self.ie_pb is not None else None],
            "ie_gb": [self.ie_gb if self.ie_gb is not None else None],
            "qh": [self.qh if self.qh is not None else None],

            "gb_energy_complex_bond": self.gb_energy_complex_bond,
            "gb_energy_complex_angle": self.gb_energy_complex_angle,
            "gb_energy_complex_dihed": self.gb_energy_complex_dihed,
            "gb_energy_complex_vdwaals": self.gb_energy_complex_vdwaals,
            "gb_energy_complex_eel": self.gb_energy_complex_eel,
            "gb_energy_complex_14vdw": self.gb_energy_complex_14vdw,
            "gb_energy_complex_14eel": self.gb_energy_complex_14eel,
            "gb_energy_complex_egb": self.gb_energy_complex_egb,
            "gb_energy_complex_esurf": self.gb_energy_complex_esurf,
            "gb_energy_complex_ggas": self.gb_energy_complex_ggas,
            "gb_energy_complex_gsolv": self.gb_energy_complex_gsolv,
            "gb_energy_complex_total": self.gb_energy_complex_total,

            "gb_energy_receptor_bond": self.gb_energy_receptor_bond,
            "gb_energy_receptor_angle": self.gb_energy_receptor_angle,
            "gb_energy_receptor_dihed": self.gb_energy_receptor_dihed,
            "gb_energy_receptor_vdwaals": self.gb_energy_receptor_vdwaals,
            "gb_energy_receptor_eel": self.gb_energy_receptor_eel,
            "gb_energy_receptor_14vdw": self.gb_energy_receptor_14vdw,
            "gb_energy_receptor_14eel": self.gb_energy_receptor_14eel,
            "gb_energy_receptor_egb": self.gb_energy_receptor_egb,
            "gb_energy_receptor_esurf": self.gb_energy_receptor_esurf,
            "gb_energy_receptor_ggas": self.gb_energy_receptor_ggas,
            "gb_energy_receptor_gsolv": self.gb_energy_receptor_gsolv,
            "gb_energy_receptor_total": self.gb_energy_receptor_total,

            "gb_energy_ligand_bond": self.gb_energy_ligand_bond,
            "gb_energy_ligand_angle": self.gb_energy_ligand_angle,
            "gb_energy_ligand_dihed": self.gb_energy_ligand_dihed,
            "gb_energy_ligand_vdwaals": self.gb_energy_ligand_vdwaals,
            "gb_energy_ligand_eel": self.gb_energy_ligand_eel,
            "gb_energy_ligand_14vdw": self.gb_energy_ligand_14vdw,
            "gb_energy_ligand_14eel": self.gb_energy_ligand_14eel,
            "gb_energy_ligand_egb": self.gb_energy_ligand_egb,
            "gb_energy_ligand_esurf": self.gb_energy_ligand_esurf,
            "gb_energy_ligand_ggas": self.gb_energy_ligand_ggas,
            "gb_energy_ligand_gsolv": self.gb_energy_ligand_gsolv,
            "gb_energy_ligand_total": self.gb_energy_ligand_total,

            "gb_energy_delta_bond": self.gb_energy_delta_bond,
            "gb_energy_delta_angle": self.gb_energy_delta_angle,
            "gb_energy_delta_dihed": self.gb_energy_delta_dihed,
            "gb_energy_delta_vdwaals": self.gb_energy_delta_vdwaals,
            "gb_energy_delta_eel": self.gb_energy_delta_eel,
            "gb_energy_delta_14vdw": self.gb_energy_delta_14vdw,
            "gb_energy_delta_14eel": self.gb_energy_delta_14eel,
            "gb_energy_delta_egb": self.gb_energy_delta_egb,
            "gb_energy_delta_esurf": self.gb_energy_delta_esurf,
            "gb_energy_delta_ggas": self.gb_energy_delta_ggas,
            "gb_energy_delta_gsolv": self.gb_energy_delta_gsolv,
            "gb_energy_delta_total": self.gb_energy_delta_total,


            "pb_energy_complex_bond": self.pb_energy_complex_bond,
            "pb_energy_complex_angle": self.pb_energy_complex_angle,
            "pb_energy_complex_dihed": self.pb_energy_complex_dihed,
            "pb_energy_complex_vdwaals": self.pb_energy_complex_vdwaals,
            "pb_energy_complex_eel": self.pb_energy_complex_eel,
            "pb_energy_complex_14vdw": self.pb_energy_complex_14vdw,
            "pb_energy_complex_14eel": self.pb_energy_complex_14eel,
            "pb_energy_complex_epb": self.pb_energy_complex_epb,
            "pb_energy_complex_enpolar": self.pb_energy_complex_enpolar,
            "pb_energy_complex_edisper": self.pb_energy_complex_edisper,
            "pb_energy_complex_ggas": self.pb_energy_complex_ggas,
            "pb_energy_complex_gsolv": self.pb_energy_complex_gsolv,
            "pb_energy_complex_total": self.pb_energy_complex_total,

            "pb_energy_receptor_bond": self.pb_energy_receptor_bond,
            "pb_energy_receptor_angle": self.pb_energy_receptor_angle,
            "pb_energy_receptor_dihed": self.pb_energy_receptor_dihed,
            "pb_energy_receptor_vdwaals": self.pb_energy_receptor_vdwaals,
            "pb_energy_receptor_eel": self.pb_energy_receptor_eel,
            "pb_energy_receptor_14vdw": self.pb_energy_receptor_14vdw,
            "pb_energy_receptor_14eel": self.pb_energy_receptor_14eel,
            "pb_energy_receptor_epb": self.pb_energy_receptor_epb,
            "pb_energy_receptor_enpolar": self.pb_energy_receptor_enpolar,
            "pb_energy_receptor_edisper": self.pb_energy_receptor_edisper,
            "pb_energy_receptor_ggas": self.pb_energy_receptor_ggas,
            "pb_energy_receptor_gsolv": self.pb_energy_receptor_gsolv,
            "pb_energy_receptor_total": self.pb_energy_receptor_total,

            "pb_energy_ligand_bond": self.pb_energy_ligand_bond,
            "pb_energy_ligand_angle": self.pb_energy_ligand_angle,
            "pb_energy_ligand_dihed": self.pb_energy_ligand_dihed,
            "pb_energy_ligand_vdwaals": self.pb_energy_ligand_vdwaals,
            "pb_energy_ligand_eel": self.pb_energy_ligand_eel,
            "pb_energy_ligand_14vdw": self.pb_energy_ligand_14vdw,
            "pb_energy_ligand_14eel": self.pb_energy_ligand_14eel,
            "pb_energy_ligand_epb": self.pb_energy_ligand_epb,
            "pb_energy_ligand_enpolar": self.pb_energy_ligand_enpolar,
            "pb_energy_ligand_edisper": self.pb_energy_ligand_edisper,
            "pb_energy_ligand_ggas": self.pb_energy_ligand_ggas,
            "pb_energy_ligand_gsolv": self.pb_energy_ligand_gsolv,
            "pb_energy_ligand_total": self.pb_energy_ligand_total,

            "pb_energy_delta_bond": self.pb_energy_delta_bond,
            "pb_energy_delta_angle": self.pb_energy_delta_angle,
            "pb_energy_delta_dihed": self.pb_energy_delta_dihed,
            "pb_energy_delta_vdwaals": self.pb_energy_delta_vdwaals,
            "pb_energy_delta_eel": self.pb_energy_delta_eel,
            "pb_energy_delta_14vdw": self.pb_energy_delta_14vdw,
            "pb_energy_delta_14eel": self.pb_energy_delta_14eel,
            "pb_energy_delta_epb": self.pb_energy_delta_epb,
            "pb_energy_delta_enpolar": self.pb_energy_delta_enpolar,
            "pb_energy_delta_edisper": self.pb_energy_delta_edisper,
            "pb_energy_delta_ggas": self.pb_energy_delta_ggas,
            "pb_energy_delta_gsolv": self.pb_energy_delta_gsolv,
            "pb_energy_delta_total": self.pb_energy_delta_total,
        }
        pd.DataFrame.from_dict(delta_g_dict).to_csv(output_file, index=False)
