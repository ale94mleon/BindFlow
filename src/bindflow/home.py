#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import bindflow
from pathlib import Path
import sys
import inspect
import platform


def home(dataDir=None, libDir=False) -> Path:
    """Return the pathname of the bindflow root directory (or a data subdirectory).
    Parameters
    ----------
    dataDir : str
        If not None, return the path to a specific data directory
    libDir : bool
        If True, return path to the lib directory
    Returns
    -------
    dir : pathlib.Path
        The directory
    Example
    -------
    .. ipython:: python

        from bindflow.home import home
        print(home())
        print(home(dataDir="gmx_ff"))
        print(home(dataDir="gmx_ff")/"amber99sb-star-ildn.ff.tar.gz")
    """

    homeDir = Path(inspect.getfile(bindflow)).parent
    try:
        if sys._MEIPASS:
            homeDir = Path(sys._MEIPASS)
    except Exception:
        pass

    if dataDir:
        return homeDir/f"data/{dataDir}"
    elif libDir:
        libdir = homeDir/f"lib{platform.system()}"
        if not libdir.exists():
            raise FileNotFoundError("Could not find libs.")
        return libdir
    else:
        return homeDir


if __name__ == "__main__":
    pass
