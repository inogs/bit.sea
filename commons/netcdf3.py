import scipy.io.netcdf as NC

def read_2d_file(filename,varname):
    ncIN = NC.netcdf_file(filename,'r')
    M2d = ncIN.variables[varname].data.copy()
    ncIN.close()
    return M2d

def write_2d_file(M2d,varname,outfile,mask,fillValue=1.e+20):
    ncOUT = NC.netcdf_file(outfile,'w')
    _, jpj, jpi= mask.shape
    ncOUT.createDimension("longitude", jpi)
    ncOUT.createDimension("latitude", jpj)
    ncvar = ncOUT.createVariable(varname, 'f', ('latitude','longitude'))
    setattr(ncvar,'fillValue'    ,fillValue)
    setattr(ncvar,'missing_value',fillValue)
    ncvar[:] = M2d
    ncOUT.close()
    
def write_3d_file(M3d,varname,outfile,mask,fillValue=1.e+20):
    '''
    mask is a Mask object
    '''
    ncOUT = NC.netcdf_file(outfile,'w')
    jpk, jpj, jpi= mask.shape
    ncOUT.createDimension("longitude", jpi)
    ncOUT.createDimension("latitude", jpj)
    ncOUT.createDimension("depth"   , jpk)
    ncvar = ncOUT.createVariable(varname, 'f', ('depth','latitude','longitude'))
    setattr(ncvar,'fillValue'    ,fillValue)
    setattr(ncvar,'missing_value',fillValue)
    ncvar[:] = M3d
    ncOUT.close()