# Copyright (c) 2015 eXact Lab srl
# Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>
import sys
import os
import numpy as np
import netCDF4
import warnings


#Mask object
from bitsea.commons.mask import Mask

class NotFoundError(Exception):
    pass

class DataExtractor(object):
    """Extracts data from a NetCDF model file.

    This class is meant to be used to read data from a model file and provide
    some facilities to access to the data
    Mask object is required.
    Nevertheless, consistency in depth between data and mask is not required.



    """
    def __init__(self, mask, filename=None, varname=None, rawdata=None, rawdatafill=np.nan, fill_value=np.nan,dimvar=3, verbose=True):
        """DataExtractor constructor.

        Args:
            - *mask*: a mask object.
            - *filename* (optional): the name of the NetCDF file containing the data.
            - *varname* (optional): the name of the variable of interest.
            - *rawdata* (optional): a Numpy array containing the data.
            - *rawdatafill* (optional): the fill value used inside rawdata
              (default: np.nan)
            - *fill_value* (optional): the value that will be used when there's
            - *dimvar* (optional) : if present, it forces the dimension of the output.
              If dimvar=2 when we read a 3D array, a 2d array of the surface will be returned.
            - * verbose * : logical, if True warnings are printed when data and mask are not consistent.

        Either rawdata or filename plus varname must be defined.
        """
        if filename is None and rawdata is None:
            raise ValueError(
                "At least one between filename and rawdata must be defined"
            )
        if filename is not None and rawdata is not None:
            raise ValueError(
                "filename and rawdata can not be submitted at the same time"
            )
        if filename is not None and varname is None:
            raise ValueError(
                "filename and varname must be submitted at the same time"
            )

        self.fill_value = fill_value

        if filename is None:
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

            # Try open the NetCDF file and search for var
            dset = netCDF4.Dataset(fn)
            if v not in dset.variables:
                raise NotFoundError("variable '%s' not found" % (varname, ))

            nc_declared_shape = np.array(dset.variables[v].shape)
            one_dimensions_ind = np.nonzero(nc_declared_shape == 1 )[0]
            if one_dimensions_ind.size == 0:  # no dimensions to remove
                self.__values = np.array(dset.variables[v])
            if one_dimensions_ind.size ==1 :
                if one_dimensions_ind[0]==0:
                    self.__values = np.array(dset.variables[v])[0,:]

            nc_declared_shape = dset.variables[v].shape

            v_slices = [slice(None) for _ in range(len(nc_declared_shape))]

            # Remove useless dimension at the beginning of the vector
            k = 0
            while k < len(v_slices):
                if nc_declared_shape[k] == 1:
                    v_slices[k] = 0
                else:
                    break
                k += 1

            new_shape = nc_declared_shape[k:]

            if len(new_shape) == 3 and dimvar == 2:
                v_slices[-3] = 0

            self.__values = np.array(dset.variables[v][tuple(v_slices)])

            self.__shape = self.__values.shape
            self.dims = len(self.__values.shape)

            for attr_name in ["missing_value","fillvalue","fillValue","FillValue"]:
                if attr_name in dset.variables[v].ncattrs():
                    fv = getattr(dset.variables[v], attr_name)
                    self.__dset_fillvalue = fv

            self.__filename = fn
            self.__varname = v

        #test dimensions
        if self.dims==3:
            if self.__shape[1:] != mask.shape[1:]:
                raise ValueError("mask must have the same shape of the data")
            data_jpk=self.__shape[0]
            mask_jpk=  mask.shape[0]
            if data_jpk > mask_jpk: # working with reduced mask
                if verbose:
                    warnings.warn('WARNING: slicing 3D field in range 0 - {}'.format(mask_jpk))
                self.__values=self.values[:mask_jpk,:,:]

            if data_jpk < mask_jpk:
                if verbose:
                    warnings.warn('WARNING: 3D file is a subset of mask domain')
                appval = self.__values
                self.__values = np.zeros(mask.shape)
                self.__values[:self.__shape[0],:,:] = appval
                self.__shape = self.__values.shape
        if self.dims==2:
            if self.__shape != mask.shape[1:]:
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

if __name__ == "__main__":
    TheMask= Mask('/g100_work/OGS_devC/Benchmark/SETUP/PREPROC/MASK/meshmask.nc')
    filename="/g100_work/OGS_prodC/OPA/V10C-prod/wrkdir/forecast/0/MODEL/AVE_FREQ_1/ave.20240229-12:00:00.N1p.nc"
    De = DataExtractor(TheMask,filename=filename,varname='N1p')
    M3d = De.values*3
    De= DataExtractor(TheMask,rawdata=M3d)

