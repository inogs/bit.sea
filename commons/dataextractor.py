# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
import sys
import os
import numpy as np
import netCDF4

#Mask object
from commons.mask import Mask

class DataExtractor(object):
    """Extracts data from a NetCDF model file.

    This class is meant to be used to read data from a model file and provide
    some facilities to access to the data
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
        self.fill_value = fill_value
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

    @property
    def values(self):
        """Get the variable values Numpy array as extracted from the file.
        """
        return self.__values

    @property
    def data_fill_value(self):
        """Get the missing data fill value as defined in the source file.
        """
        return self.__dset_fillvalue

    @property
    def filled_values(self):
        """Get a copy of the values but fill the missing values with
        self.fill_value.
        """
        output = np.copy(self.__values)
        output[self.__values == self.__dset_fillvalue] = self.fill_value
        return output