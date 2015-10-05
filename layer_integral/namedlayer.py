import sys
import os
import numpy as np
import netCDF4 as nc

#Kludge added for testing purpose TO BE REMOVED
sys.path.append("/home/skyglobe/git-repos/bit.sea/matchup")

#Layer Object
from matchup.matchup import Layer

class NamedLayer(Layer):
    """
    Holds layer information taken from a NetCDF file
    """
    def __init__(self, top, bottom, var, ncfilename):
        t = float(top)
        b = float(bottom)
        if t < b:
            self.__top = t
            self.__bottom = b
        else:
            raise ValueError("top must be above of bottom")
        self._var = str(var)
        self._ncfilename = str(ncfilename)

    @property
    def top(self):
        return self.__top

    @property
    def bottom(self):
        return self.__bottom

