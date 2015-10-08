# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
import sys
import os
import numpy as np
import netCDF4

#Layer Object
from layer import Layer

#Mask object
from mask import Mask

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
                raise ValueError("variable '%s' not found" % (var, ))
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

    def get_integral(self, timestep=0, fill_value=np.nan):
        """Returns a 2D NumPy array with the integral over depth.
        """
        #Find Z indices
        top_index = np.where(self._mask.zlevels >= self.__top)[0][0]
        bottom_index = np.where(self._mask.zlevels < self.__bottom)[0][-1]
        if top_index == bottom_index:
            #Just one layer so we return the sliced data
            output = np.array(dset.variables[self.__varname][timestep,top_index,:,:])
            fv = dset.variables[self.__varname].missing_value
            output[output == fv] = fill_value
            return output
        #Workaround for Python ranges
        bottom_index += 1
        #Build local mask matrix
        lmask = np.array(self._mask.mask[top_index:bottom_index,:,:], dtype=np.double)
        #Build dz matrix
        dzm = np.ones_like(lmask, dtype=np.double)
        j = 0
        for i in range(top_index, bottom_index):
            dzm[j,:,:] = self._mask.dz[i]
            j += 1
        #Get the slice of the values
        v = np.array(dset.variables[self.__varname][timestep,top_index:bottom_index,:,:])
        #Build integral matrix (2D)
        integral = (v * dzm * lmask).sum(axis=0)
        #Build height matrix (2D)
        height = (dzm * lmask).sum(axis=0)
        indexmask = [height > 0]
        #Build output matrix (2D)
        output = np.full_like(integral, fill_value, dtype=np.double)
        output[indexmask] = integral[indexmask] / height[indexmask]
        return output
