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

class DataExtractor(object):
    """Extracts data from a NetCDF model file.

    This class is meant to be used to read data from a model file and provide
    some facilities to access to the post-processed data
    """
    def __init__(self, varname, filename, mask, fill_value=np.nan):
        """DataExtractor constructor.

        Args:
            - *varname*: the name of the variable of interest.
            - *filename*: the name of the NetCDF file containing the data.
            - *mask*: a mask object.
            _ *fill_value* (optional): the value that will be used when there's
              no data (default: np.nan)
        """
        fn = str(filename)
        v = str(varname)
        self.__fill_value = fill_value
        try:
            #Try open the NetCDF file and search for var
            dset = netCDF4.Dataset(fn)
            if not v in dset.variables:
                raise ValueError("variable '%s' not found" % (varname, ))
            else:
                self.__shape = dset.variables[v].shape
                self.__values = np.array(dset.variables[v])
                fv = dset.variables[v].missing_value
                self.__dset_fillvalue = fv
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

    def get_layer_average(self, layer, timestep=0):
        """Returns a 2D NumPy array with the average weighted over depth.

        Args:
            - *layer*: a Layer object that contains the vertical boundaries for
              the computation.
            - *timestep* (optional): the index of the first dimension (time) in
              the model data. Default: 0.
        """
        if not isinstance(layer, (Layer,)):
            raise ValueError("layer must be a Layer object")
        #Find Z indices
        top_index = np.where(self._mask.zlevels >= layer.top)[0][0]
        bottom_index = np.where(self._mask.zlevels < layer.bottom)[0][-1]
        if top_index == bottom_index:
            #Just one layer so we return the sliced data
            output = np.array(self.__values[timestep,top_index,:,:])
            output[output == self.__dset_fillvalue] == self.__fill_value
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
        v = np.array(self.__values[timestep,top_index:bottom_index,:,:])
        #Build integral matrix (2D)
        integral = (v * dzm * lmask).sum(axis=0)
        #Build height matrix (2D)
        height = (dzm * lmask).sum(axis=0)
        indexmask = [height > 0]
        #Build output matrix (2D)
        output = np.full_like(integral, self.__fill_value, dtype=np.double)
        output[indexmask] = integral[indexmask] / height[indexmask]
        return output
