import netCDF4 as NC
import os
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


def readfile(filename, var):
    '''
    Generic file reader
    '''
    D = NC.Dataset(filename,"r")
    VAR = np.array(D. variables[var])
    D.close()
    return VAR
    

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

#     if os.path.exists(outfile):
#         ncOUT=NC.Dataset(outfile,'a')
#         print "appending ", varname, " in ", outfile
#     else:
    ncOUT = NC.Dataset(outfile,'w')
    jpk, jpj, jpi= mask.shape
    ncOUT.createDimension("longitude", jpi)
    ncOUT.createDimension("latitude", jpj)
    ncOUT.createDimension("depth"   , jpk)
    
    dims = (depth_dimension_name(ncOUT),lat_dimension_name(ncOUT),lon_dimension_name(ncOUT))
    ncvar = ncOUT.createVariable(varname, 'f', dims, zlib=True, fill_value=1.0e+20) 
    setattr(ncvar,'fillValue'    ,fillValue)
    setattr(ncvar,'missing_value',fillValue)
    ncvar[:] = M3d
    ncOUT.close()