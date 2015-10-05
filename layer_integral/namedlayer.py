import sys
import os
import numpy as np
import netCDF4

#Kludge added for testing purpose TO BE REMOVED
sys.path.append("/home/skyglobe/git-repos/bit.sea/matchup")

#Layer Object
from matchup.matchup import Layer

class NamedLayer(Layer):
    """
    Holds layer information taken from a NetCDF file
    """
    def __init__(self, top, bottom, varname, filename):
        t = float(top)
        b = float(bottom)
        if t < b:
            self.__top = t
            self.__bottom = b
        else:
            raise ValueError("top must be above of bottom")

        fn = str(filename)
        v = str(varname)
        try:
            #Try open the NetCDF file and search for var
            dset = netCDF4.Dataset(fn)
            if not v in dset.variables:
                raise ValueError("variable '%s' not found" % (var, ))
            dset.close()
            self.__filename = fn
            self.__varname = v
        except:
            raise


    @property
    def top(self):
        return self.__top

    @property
    def bottom(self):
        return self.__bottom

    @property
    def filename(self):
        return self.__filename

    @property
    def varname(self):
        return self.__varname
