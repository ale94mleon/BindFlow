#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import bindflow
import os
import sys
import inspect
import platform


def home(dataDir=None, libDir=False):
    """Return the pathname of the bindflow root directory (or a data subdirectory).
    Parameters
    ----------
    dataDir : str
        If not None, return the path to a specific data directory
    libDir : bool
        If True, return path to the lib directory
    Returns
    -------
    dir : str
        The directory
    Example
    -------
    .. ipython:: python

        from bindflow.home import home
        import os
        print(home())
        print(home(dataDir="gmx_ff"))
        os.path.join(home(dataDir="gmx_ff"),"amber99sb-star-ildn.ff.tar.gz")
    """

    homeDir = os.path.dirname(inspect.getfile(bindflow))
    try:
        if sys._MEIPASS:
            homeDir = sys._MEIPASS
    except Exception:
        pass

    if dataDir:
        return os.path.join(homeDir, "data", dataDir)
    elif libDir:
        libdir = os.path.join(homeDir, "lib", platform.system())
        if not os.path.exists(libdir):
            raise FileNotFoundError("Could not find libs.")
        return libdir
    else:
        return homeDir


if __name__ == "__main__":
    pass
