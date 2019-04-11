from __future__ import print_function
import scipy.io.netcdf as NC
import os

def read_2d_file(filename,varname):
    '''
    Arguments:
       * filename * : string (an existing path file)
       * varname  * : string (an existing variable in path file)
    
    Hypothesis: in filename there is a variable called 'varname'
    having dimensions (time,lat,lon) or (lat,lon)
    Returns : 
        3D numpy array
    '''
    ncIN = NC.netcdf_file(filename,'r')
    nDim = len(ncIN.variables[varname].dimensions)
    if nDim==3 : VAR   = ncIN.variables[varname].data[0,:,:].copy()
    if nDim==2 : VAR   = ncIN.variables[varname].data.copy()
    ncIN.close()
    return VAR


def read_3d_file(filename, varname):
    '''
    Arguments:
       * filename * : string (an existing path file)
       * varname  * : string (an existing variable in path file)
    
    Hypothesis: in filename there is a variable called 'varname'
    having dimensions (time,depth,lat,lon) or (depth,lat,lon)
    Returns : 
        3D numpy array
    '''
    ncIN = NC.netcdf_file(filename,'r')
    nDim = len(ncIN.variables[varname].dimensions)
    if nDim==4 : VAR   = ncIN.variables[varname].data[0,:,:,:].copy()
    if nDim==3 : VAR   = ncIN.variables[varname].data.copy()
    ncIN.close()
    return VAR
    
def lon_dimension_name(ncObj):
    '''
    Argument:
        ncObj : a NetCDF object, got by NC.netcdf_file()
    '''
    for dimname in ['lon','longitude']:
        if ncObj.dimensions.has_key(dimname):
            break
    return dimname

def lat_dimension_name(ncObj):
    '''
    Argument:
        ncObj : a NetCDF object, got by NC.netcdf_file()
    '''
    for dimname in ['lat','latitude']:
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

    
def write_2d_file(M2d,varname,outfile,mask,fillValue=1.e+20):
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
        ncOUT=NC.netcdf_file(outfile,'a')
        print("appending ", varname, " in ", outfile)
    else:
        ncOUT = NC.netcdf_file(outfile,'w')
        jpk, jpj, jpi= mask.shape
        ncOUT.createDimension("longitude", jpi)
        ncOUT.createDimension("latitude" , jpj)
        ncOUT.createDimension("depth"    , jpk)

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
        
        
    ncvar = ncOUT.createVariable(varname, 'f', (lat_dimension_name(ncOUT),lon_dimension_name(ncOUT)))
    setattr(ncvar,'fillValue'    ,fillValue)
    setattr(ncvar,'missing_value',fillValue)
    setattr(ncvar,'coordinates'  ,'latitude longitude')
    ncvar[:] = M2d
    setattr(ncOUT,'latitude_min', 30.0)
    setattr(ncOUT,'latitude_max', 46.0)
    setattr(ncOUT,'longitude_min',-6.0)
    setattr(ncOUT,'longitude_max',37.0)
    ncOUT.close()
    
def write_3d_file(M3d,varname,outfile,mask,fillValue=1.e+20):
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
        ncOUT=NC.netcdf_file(outfile,'a')
        print("appending ", varname, " in ", outfile)
    else:
        ncOUT = NC.netcdf_file(outfile,'w')
        jpk, jpj, jpi= mask.shape
        ncOUT.createDimension("longitude", jpi)
        ncOUT.createDimension("latitude", jpj)
        ncOUT.createDimension("depth"   , jpk)
    

    ncvar = ncOUT.createVariable(varname, 'f', (depth_dimension_name(ncOUT),lat_dimension_name(ncOUT),lon_dimension_name(ncOUT)))
    setattr(ncvar,'fillValue'    ,fillValue)
    setattr(ncvar,'missing_value',fillValue)
    ncvar[:] = M3d
    ncOUT.close()

