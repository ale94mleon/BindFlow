import copy
import os
from typing import Union

import importlib_resources
from GMXMMPBSA.input_parser import InputFile
from GMXMMPBSA.input_parser import input_file as _InternalInputFile

PathLike = Union[os.PathLike, str, bytes]
_bindflow_template_resources = importlib_resources.files("bindflow")


class GMXMMPBSAInputMaker:
    """Class to handle gmx_MMPBSA input files

    Example
    -------

    .. ipython:: python

        import sys
        input_file = GMXMMPBSAInputMaker(pb={"ipb": 123}, general={"sys_name": "abc"})
        # In real example, just provided the file path
        input_file.write(sys.stdout)

    """
    def __init__(self, **kwargs):
        """Here we check the type of calculation to perform based on the provided kwargs.
        If no kwargs was provided a default Poison-Boltzmann (pb) calculation is performed
        based on the internal BindFlow templates. If you would like to perform a General Boltzmann (gb)
        calculation you must provided the keyword gb={}. You can specify parameters inside the dictionary
        or leave it empty and the internal default BindFlow parameters will be used. See that if the values
        of each passed keyword is not a dictionary, a ValueError will be raised. E.g. gb=True is an invalid
        keyword.
        """
        for key, value in kwargs.items():
            if not isinstance(value, dict):
                raise ValueError(f"{key} = {value} is an invalid keyword. The value must be a dictionary")
        self.__mmpbsa_in_opts = kwargs
        self.__do_gb = True if "gb" in self.__mmpbsa_in_opts.keys() else False
        # Do pb by default
        self.__do_pb = True if "pb" in self.__mmpbsa_in_opts.keys() or not self.__do_gb else False
        self.__template = self.__load_from_template()
        self.__update_user_specific_fields()

    def __load_from_template(self) -> InputFile:
        """Loading internal BindFLow templates on demand.

        Returns
        -------
        InputFile
            An instance of `GMXMMPBSA.input_parser.InputFile` with
            the namelist initiated internally by  GMXMMPBSA.
        """
        input_copy = copy.deepcopy(_InternalInputFile)
        if self.__do_pb and self.__do_gb:
            template_path = _bindflow_template_resources.joinpath("mmpbsa_in/templates", "pb_gb.in")
            input_copy.Parse(template_path)
            return input_copy
        elif self.__do_pb:
            template_path = _bindflow_template_resources.joinpath("mmpbsa_in/templates", "pb.in")
            input_copy.Parse(template_path)
            return input_copy
        elif self.__do_gb:
            template_path = _bindflow_template_resources.joinpath("mmpbsa_in/templates", "gb.in")
            input_copy.Parse(template_path)
            return input_copy

    def __update_user_specific_fields(self) -> None:
        """Helper function to pass user defined parameters
        passed during initialization of the class

        Raises
        ------
        ValueError
            In case of invalid parameter
        """
        if self.__mmpbsa_in_opts:
            for key in self.__mmpbsa_in_opts.keys():
                if self.__mmpbsa_in_opts[key] is None or self.__mmpbsa_in_opts[key] is True:
                    continue
                for parameter in self.__mmpbsa_in_opts[key].keys():
                    if parameter in self.__template.namelists[key].variables.keys():
                        self.__template.namelists[key].variables[parameter].SetValue(self.__mmpbsa_in_opts[key][parameter])
                    else:
                        raise ValueError(f"The parameter {key}/{parameter} for the MMPBSA/MMGBSA calculation is unknown. "
                                         "Check this list https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/input_file/ "
                                         "for possible input options.")

    def write(self, out_path: PathLike):
        outputs = ["general"]
        if self.__do_gb:
            outputs.append("gb")
        if self.__do_pb:
            outputs.append("pb")
        self.__template.print_contents(out_path, outputs)


if __name__ == "__main__":
    pass
