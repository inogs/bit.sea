from postproc import masks
import scipy.io.netcdf as NC
import numpy as np



NativeMesh = masks.SatOrigMesh
fillValue = -999.0

    
def readfromfile(filename):
    '''
    returns CHL
    ''' 
    ncIN = NC.netcdf_file(filename,'r')
    CHL_IN = ncIN.variables['CHL'].data[0,:,:].copy()
    ncIN.close()
    return CHL_IN

def readClimatology(filename):
    ncIN = NC.netcdf_file(filename,'r') 
    MEAN = ncIN.variables['Mean'].data #(366, 253, 733)
    STD  = ncIN.variables['Mean'].data 
    ncIN.close()
    return MEAN,STD

def dumpfile(self,filename, CHL):
    ncOUT  = NC.netcdf_file(filename,'w')
    ncOUT.createDimension('time', 1)
    ncOUT.createDimension('depth',1)
    ncOUT.createDimension('lon',NativeMesh.jpi)
    ncOUT.createDimension('lat',NativeMesh.jpj)
    
    ncvar = ncOUT.createVariable('CHL', 'f', ('lat','lon'))
    ncvar[:] = CHL
    setattr(ncvar, 'missing_value', fillValue)
    setattr(ncvar, '_FillValue', fillValue)
    ncvar = ncOUT.createVariable('depth','f',('depth',))
    ncvar[:] = 1.47210180759
    ncvar = ncOUT.createVariable('time','f',('time',))
    ncvar[:] = 1.0
    ncvar = ncOUT.createVariable('lat','f',('lat',))
    ncvar[:] = NativeMesh.lat
    ncvar = ncOUT.createVariable('lon','f',('lon',))
    ncvar[:] = NativeMesh.lon
    ncOUT.close()   
def dumpV4file(self,filename,M):
    a=0
def dumpV1file(self,filename,M):
    a=0
    
def logAverager(self,M):
    '''
    Inner matrix M has dimensions (nFrames, jpj, jpi )
    Performs average passing through natural logarithm
    Mean  = exp(ln(values).mean() ) 
     
    '''
    nFrames,jpj,jpi = M.shape
    assert jpj == NativeMesh.jpj
    assert jpi == NativeMesh.jpi
    CHL_OUT = np.ones((NativeMesh.jpj,NativeMesh.jpi),np.float32) * fillValue
    
    for i in range(NativeMesh.jpi):
        for j in range(NativeMesh.jpj):
            l = M[:,j,i]
            goodValues = l != fillValue
            if np.any(goodValues):
                count = goodValues.sum()
                LOGmean = np.log(l[goodValues]).sum() / count
                CHL_OUT[j,i] = np.exp(LOGmean)
    return CHL_OUT
    
    