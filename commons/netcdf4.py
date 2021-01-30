from __future__ import print_function
import netCDF4 as NC
import os
import numpy as np
def lon_dimension_name(ncObj):
    '''
    Argument:
        ncObj : a NetCDF object, got by NC.netcdf_file()
    '''
    for dimname in ['lon','longitude','x']:
        if ncObj.dimensions.has_key(dimname):
            break
    return dimname

def lat_dimension_name(ncObj):
    '''
    Argument:
        ncObj : a NetCDF object, got by NC.netcdf_file()
    '''
    for dimname in ['lat','latitude','y']:
        if ncObj.dimensions.has_key(dimname):
            break
    return dimname

def depth_dimension_name(ncObj):
    '''
    Argument:
        ncObj : a NetCDF object, got by NC.netcdf_file()
    '''
    for dimname in ['depth','z']:
        if ncObj.dimensions.has_key(dimname):
            break
    return dimname


def readfile(filename, var):
    '''
    Generic file reader
    '''
    D = NC.Dataset(filename,"r")
    VAR = np.array(D. variables[var])
    D.close()
    return VAR
    

def write_3d_file(M3d,varname,outfile,mask,fillValue=1.e+20, compression=False):
    '''
    Dumps a 3D array in a NetCDF file.


    Arguments:
    * M3d       * the 3D array to dump
    * varname   * the variable name on NetCDF file
    * outfile   * file that will be created. If it is an existing file,
                  it will be opened in 'append' mode.
    * mask      * a mask object consistent with M3d array
    * fillvalue * (optional) value to set missing_value attribute.

    When the file is opened in 'append' mode this method tries to adapt to
    existing dimension names (for example it works both with 'lon' or 'longitude')

    Does not return anything.
    '''

    if os.path.exists(outfile):
        ncOUT=NC.Dataset(outfile,'a')
        print("appending ", varname, " in ", outfile)
        variable_exist= ncOUT.variables.has_key(varname)
        if variable_exist:
            ncvar=ncOUT.variables[varname]
        else:
            dims = (depth_dimension_name(ncOUT),lat_dimension_name(ncOUT),lon_dimension_name(ncOUT))
            ncvar = ncOUT.createVariable(varname, 'f', dims, zlib=compression, fill_value=fillValue)
            setattr(ncvar,'fillValue'    ,fillValue)
            #setattr(ncvar,'missing_value',fillValue)
    else:
        ncOUT = NC.Dataset(outfile,'w')
        
        jpk, jpj, jpi= mask.shape
        ncOUT.createDimension("longitude", jpi)
        ncOUT.createDimension("latitude", jpj)
        ncOUT.createDimension("depth"   , jpk)
        dims = (depth_dimension_name(ncOUT),lat_dimension_name(ncOUT),lon_dimension_name(ncOUT))
        ncvar = ncOUT.createVariable(varname, 'f', dims, zlib=compression, fill_value=fillValue)
        setattr(ncvar,'fillValue'    ,fillValue)

    ncvar[:] = M3d
    ncOUT.close()

def write_2d_file(M2d,varname,outfile,mask,fillValue=1.e+20, compression=False):
    '''
    Dumps a 2D array in a NetCDF file.


    Arguments:
    * M2d       * the 2D array to dump
    * varname   * the variable name on NetCDF file
    * outfile   * file that will be created. If it is an existing file,
                  it will be opened in 'append' mode.
    * mask      * a mask object consistent with M2d array
    * fillvalue * (optional) value to set missing_value attribute.

    When the file is opened in 'append' mode this method tries to adapt to
    existing dimension names (for example it works both with 'lon' or 'longitude')

    Does not return anything.'''

    if os.path.exists(outfile):
        ncOUT=NC.Dataset(outfile,'a')
        print("appending ", varname, " in ", outfile)
        variable_exist= ncOUT.variables.has_key(varname)
        if variable_exist:
            ncvar2d=ncOUT.variables[varname]
        else:
            dims = (lat_dimension_name(ncOUT),lon_dimension_name(ncOUT))
            ncvar2d = ncOUT.createVariable(varname, 'f', dims, zlib=compression, fill_value=fillValue)
            setattr(ncvar2d,'fillValue'    ,fillValue)

    else:
        ncOUT = NC.Dataset(outfile,'w')
        _, jpj, jpi= mask.shape
        ncOUT.createDimension("longitude", jpi)
        ncOUT.createDimension("latitude" , jpj)

        ncvar = ncOUT.createVariable('longitude','f', ('longitude',))
        setattr(ncvar, 'units'        ,'degrees_east')
        setattr(ncvar,'long_name'    ,'longitude')
        setattr(ncvar, 'standard_name','longitude')
        setattr(ncvar, 'axis'         ,'X')
        setattr(ncvar, 'valid_min'    , -5.5625)
        setattr(ncvar, 'valid_max'    , 36.25)
        setattr(ncvar, '_CoordinateAxisType',"Lon" )
        ncvar[:] = mask.xlevels[0,:]

        ncvar = ncOUT.createVariable( 'latitude','f', ('latitude',))
        setattr(ncvar, 'units'        ,'degrees_north')
        setattr(ncvar,'long_name'    ,'latitude')
        setattr(ncvar,'standard_name','latitude')
        setattr(ncvar, 'axis'         ,'Y')
        setattr(ncvar,'valid_min'    , 30.1875)
        setattr(ncvar, 'valid_max'    , 45.9375)
        setattr(ncvar, '_CoordinateAxisType',"Lat" )
        ncvar[:] = mask.ylevels[:,0]

        dims = (lat_dimension_name(ncOUT),lon_dimension_name(ncOUT))
        ncvar2d = ncOUT.createVariable(varname, 'f', dims, zlib=compression, fill_value=fillValue)
        setattr(ncvar2d,'fillValue'    ,fillValue)
        setattr(ncvar2d,'coordinates'  ,'latitude longitude')
    ncvar2d[:] = M2d
    setattr(ncOUT,'latitude_min', 30.0)
    setattr(ncOUT,'latitude_max', 46.0)
    setattr(ncOUT,'longitude_min',-6.0)
    setattr(ncOUT,'longitude_max',37.0)
    ncOUT.close()


def dimfile(filename, varname):
    dset = NC.Dataset(filename)
    var_obj=dset.variables[varname]
    ndims=len(var_obj.dimensions)
    truedims = ndims
    if 'time' in var_obj.dimensions:
        truedims =ndims-1
    if 'time_counter' in var_obj.dimensions:
        truedims =ndims-1     
    dset.close()
    return truedims
