import GMXMMPBSA.input_parser as inp
import copy
import importlib_resources


class gmxmmpbsa_input_file:
    __bindflow_template_resources = importlib_resources.files("bindflow")
    def __init__(self, mmpbsa_in_opts: dict = {"pb":None}):
        self.__mmpbsa_in_opts = mmpbsa_in_opts
        self.__do_gb = True if "gb" in mmpbsa_in_opts.keys() else False
        self.__do_pb = True if "pb" in mmpbsa_in_opts.keys() else False
        self.__template = self.__load_from_template()
        self.__update_user_specific_fields()


    def __create_general_template():
        return copy.deepcopy(inp.input_file)


    def __load_from_template(self):
        input_copy = gmxmmpbsa_input_file.__create_general_template()
        if self.__do_pb and self.__do_gb:
            template_path = gmxmmpbsa_input_file.__bindflow_template_resources.joinpath("mmpbsa_in/templates", "pb_gb.in")
            input_copy.Parse(template_path)
            return input_copy
        elif self.__do_pb:
            template_path = gmxmmpbsa_input_file.__bindflow_template_resources.joinpath("mmpbsa_in/templates", "pb.in")
            input_copy.Parse(template_path)
            return input_copy
        elif self.__do_gb:
            template_path = gmxmmpbsa_input_file.__bindflow_template_resources.joinpath("mmpbsa_in/templates", "gb.in")
            input_copy.Parse(template_path)
            return input_copy


    def __update_user_specific_fields(self):
        if self.__mmpbsa_in_opts is None:
            return
        for key in self.__mmpbsa_in_opts.keys():
            if self.__mmpbsa_in_opts[key] is None or self.__mmpbsa_in_opts[key] == True:
                continue
            for parameter in self.__mmpbsa_in_opts[key].keys():
                if parameter in self.__template.namelists[key].variables.keys():
                    self.__template.namelists[key].variables[parameter].SetValue(self.__mmpbsa_in_opts[key][parameter])
                else:
                    raise ValueError(f"The parameter {key}/{parameter} for the MMPBSA/MMGBSA calculation is unknown. "
                                     "Check this list https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/input_file/ "
                                     "for possible input options.")


    def store_template(self, out_path):
        outputs = ["general"]
        if self.__do_gb:
            outputs.append("gb")
        if self.__do_pb:
            outputs.append("pb")
        self.__template.print_contents(out_path, outputs)



if __name__ == "__main__":
    import sys
    a = gmxmmpbsa_input_file({"pb":{"ipb": 123}, "general": {"sys_name": "abc"}})
    a.store_template(sys.stdout)