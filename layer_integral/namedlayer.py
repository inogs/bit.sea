# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
import sys
import os
import numpy as np
import netCDF4
import matplotlib.pyplot as pl

#Layer Object
from layer import Layer

#Mask object
from commons.mask import Mask

#Data extractor
from commons.dataextractor import DataExtractor

class NamedLayer(Layer):
    """
    Holds layer information taken from a NetCDF file
    """
    def __init__(self, top, bottom, varname, filename, mask):
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
                raise ValueError("variable '%s' not found" % (varname, ))
            else:
                self.__shape = dset.variables[v].shape
            dset.close()
            self.__filename = fn
            self.__varname = v
        except:
            raise

        if not isinstance(mask, (Mask,)):
            raise ValueError("mask must be a Mask object")
        else:
            #test dimensions
            if self.__shape != mask.shape:
                if self.__shape[1:] != mask.shape:
                    raise ValueError("mask must have the same shape of the data")
        #Preserve mask reference
        self._mask = mask

    def __str__(self):
        return "Layer %d-%d m (%s)" % (self.__top, self.__bottom, self.__varname)

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

    @property
    def shape(self):
        return self.__shape

    def get_raw_profile_data(self, x, y, timestep=0, fill_value=np.nan):
        indices = [np.logical_and(self._mask.zlevels >= self.__top, self._mask.zlevels < self.__bottom)]
        dset = netCDF4.Dataset(self.__filename)
        v = dset.variables[self.__varname][timestep,:,y,x]
        fill_val = dset.variables[self.__varname].missing_value
        dset.close()
        v = np.array(v[indices])
        v[v == fill_val] = fill_value
        return np.array(self._mask.zlevels[indices]), v

    def plot(self, timestep=0, fill_value=np.nan):
        de = DataExtractor(self.varname, self.filename, self._mask, fill_value)
        data = de.get_layer_average(self, timestep)
        plot_data = np.flipud(data)
        pl.imshow(plot_data)
        pl.colorbar()
        pl.show()
