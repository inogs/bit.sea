# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
import sys
import os
import numpy as np
import netCDF4

#Mask object
from commons.mask import Mask

class NotFoundError(Exception):
    pass

class DataExtractor(object):
    """Extracts data from a NetCDF model file.

    This class is meant to be used to read data from a model file and provide
    some facilities to access to the data
    """
    def __init__(self, mask, filename=None, varname=None, rawdata=None, rawdatafill=np.nan, fill_value=np.nan,dimvar=3):
        """DataExtractor constructor.

        Args:
            - *mask*: a mask object.
            - *filename* (optional): the name of the NetCDF file containing the data.
            - *varname* (optional): the name of the variable of interest.
            - *rawdata* (optional): a Numpy array containing the data.
            - *rawdatafill* (optional): the fill value used inside rawdata
              (default: np.nan)
            _ *fill_value* (optional): the value that will be used when there's

        Either rawdata or filename plus varname must be defined.
        """
        if not ((filename is None) ^ (rawdata is None)):
            raise ValueError("Either rawdata or filename plus varname must be defined")
        elif (not (filename is None)) and (varname is None):
            raise ValueError("filename and varname must be both defined")

        self.fill_value = fill_value

        if (filename is None):
            #Use rawdata
            self.__shape = rawdata.shape
            self.dims = len(rawdata.shape)
            self.__values = rawdata
            self.__dset_fillvalue = rawdatafill
            self.__filename = 'rawdata'
            self.__varname = 'rawdata'
        else:
            fn = str(filename)
            v = str(varname)
            try:
                #Try open the NetCDF file and search for var
                dset = netCDF4.Dataset(fn)
                if not v in dset.variables:
                    raise NotFoundError("variable '%s' not found" % (varname, ))
                else:
                    nc_declared_shape = np.array(dset.variables[v].shape)
                    one_dimensions_ind = np.nonzero(nc_declared_shape == 1 )[0]
                    if one_dimensions_ind.size == 0: # no dimensions to remove
                        self.__values = np.array(dset.variables[v])
                    if one_dimensions_ind.size ==1 :
                        if one_dimensions_ind[0]==0:
                            self.__values = np.array(dset.variables[v])[0,:]

                    if one_dimensions_ind.size==2 :
                        if one_dimensions_ind[0] ==0 :
                            if one_dimensions_ind[1] ==0 :
                                self.__values = np.array(dset.variables[v])[0,0,:]

                    self.__shape = self.__values.shape
                    self.dims  = len(self.__values.shape)
                    
                    for attr_name in ["missing_value","fillvalue","fillValue","FillValue"]:
                        if attr_name in dset.variables[v].ncattrs():
                            fv = getattr(dset.variables[v], attr_name)
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
            if self.dims==3:
                if self.__shape[1:] != mask.shape[1:]: raise ValueError("mask must have the same shape of the data")
                if self.__shape[0] > mask.shape[0]: raise ValueError("mask must have the same shape of the data")
                if self.__shape[0] < mask.shape[0]:
                    print('WARNING: loading field limited to a certain level NOT complete full 3d dataset')
                    appval = self.__values
                    self.__values = np.zeros(mask.shape)
                    self.__values[:self.__shape[0],:,:] = appval
                    self.__shape = self.__values.shape
            if self.dims==2:
                if self.__shape != mask.shape[1:] : raise ValueError("mask must have the same shape of the data")
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
