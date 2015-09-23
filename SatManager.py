from postproc import masks
import scipy.io.netcdf as NC
import numpy as np



NativeMesh = masks.SatOrigMesh
V4 = masks.V4mesh
fillValue = -999.0

    
def readfromfile(filename):
    '''
    returns CHL
    ''' 
    ncIN = NC.netcdf_file(filename,'r')
    varObj = ncIN.variables['CHL']
    ndims = len(varObj.shape)
    if ndims==2: CHL_IN=varObj.data.copy()
    if ndims==3: CHL_IN=varObj.data[0,:,:].copy()    
    ncIN.close()
    return CHL_IN

def readClimatology(filename):
    ncIN = NC.netcdf_file(filename,'r') 
    MEAN = ncIN.variables['Mean'].data #(366, 253, 733)
    STD  = ncIN.variables['Std'].data 
    ncIN.close()
    return MEAN,STD

def dumpfile(filename, CHL):
    ncOUT  = NC.netcdf_file(filename,'w')
    ncOUT.createDimension('time', 1)
    ncOUT.createDimension('depth',1)
    ncOUT.createDimension('lon',NativeMesh.jpi)
    ncOUT.createDimension('lat',NativeMesh.jpj)
    
    ncvar = ncOUT.createVariable('CHL', 'f', ('time','lat','lon'))
    chl = np.zeros((1,NativeMesh.jpj,NativeMesh.jpi),np.float32)
    chl[0,:,:] = CHL
    ncvar[:] = chl
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
    
def dumpV4file(filename,CHL):
    ncOUT  = NC.netcdf_file(filename,'w')
    ncOUT.createDimension('time', 1)
    ncOUT.createDimension('depth',1)
    ncOUT.createDimension('lon',V4.jpi)
    ncOUT.createDimension('lat',V4.jpj)
    
    ncvar = ncOUT.createVariable('CHL', 'f', ('lat','lon'))
    chl = np.ones(V4.jpj, V4.jpi) * fillValue
    chl[1:,:] = CHL[0:-1:,11:]
    ncvar[:] = chl
    setattr(ncvar, 'missing_value', fillValue)
    setattr(ncvar, '_FillValue', fillValue)
    ncvar = ncOUT.createVariable('depth','f',('depth',))
    ncvar[:] = 1.47210180759
    ncvar = ncOUT.createVariable('time','f',('time',))
    ncvar[:] = 1.0
    ncvar = ncOUT.createVariable('lat','f',('lat',))
    ncvar[:] = V4.lat
    ncvar = ncOUT.createVariable('lon','f',('lon',))
    ncvar[:] = V4.lon
    ncOUT.close()
    

def interpOnV1(CHL_16,bathy_threshold = 0):
    '''
    Performs interpolation on V1 grid using the following algorithm
    For a given pixel of V1 mesh, we select the four pixels of the OrigSat mesh 
    corresponding exactly to it.
    Then, we look at the bathymetry in these pixels, and filter pixels having bathymetry greater than
    bathy_threshold. 
    Returned value is the average of the filtered pixels.
    Condition to return values are:
     - at least 3 pixels having bathymetry < bathy_threshold
     - 2 pixels if they are on a diagonal
    
    '''
    V1 = masks.V1mesh
    Bathyfile="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/MODIS/Bathy_MODIS.nc"
    ncIN = NC.netcdf_file(Bathyfile,'r')
    BATHY=ncIN.variables['bathy'].data
    ncIN.close()
    
    PREVLON = np.array( [getnextIndex(NativeMesh.lon, V1.lon[i]) for i in range(V1.jpi)], np.int32) -1
    PREVLAT = np.array( [getnextIndex(NativeMesh.lat, V1.lat[i]) for i in range(V1.jpj)], np.int32) -1

    CHL_8 = np.ones((V1.jpj,V1.jpi),np.float32) * fillValue
    
    for i in range(V1.jpi):
        ji = PREVLON[i]
        for j in range(V1.jpj):
            jj = PREVLAT[j] 
            
            riq =CHL_16[jj:jj+2,ji:ji+2]
            bat = BATHY[jj:jj+2,ji:ji+2]
            
            batgoods = bat < bathy_threshold
            condition = False
            
            if batgoods.sum() > 3:
                condition = True               
            if batgoods.sum() == 2:
                diag1 = batgoods.diagonal()
                diag2 = [batgoods[0,1], batgoods[1,0]]
                if (diag1.all()) | (diag2.all()) :
                    condition=True
           
            if condition:
                selected = riq[batgoods]
                CHL_8[j,i] = selected[selected>fillValue].mean()
    return CHL_8    
       
def dumpV1file(filename,CHL):
 
    V1 = masks.V1mesh
    ncOUT  = NC.netcdf_file(filename,'w')
    ncOUT.createDimension('time', 1)
    ncOUT.createDimension('depth',1)
    ncOUT.createDimension('lon',V1.jpi)
    ncOUT.createDimension('lat',V1.jpj)
    
    ncvar = ncOUT.createVariable('CHL', 'f', ('lat','lon'))
    ncvar[:] = CHL
    setattr(ncvar, 'missing_value', fillValue)
    setattr(ncvar, '_FillValue', fillValue)
    ncvar = ncOUT.createVariable('depth','f',('depth',))
    ncvar[:] = 1.47210180759
    ncvar = ncOUT.createVariable('time','f',('time',))
    ncvar[:] = 1.0
    ncvar = ncOUT.createVariable('lat','f',('lat',))
    ncvar[:] = V1.lat
    ncvar = ncOUT.createVariable('lon','f',('lon',))
    ncvar[:] = V1.lon
    ncOUT.close()                  

    
def logAverager(M):
    '''
    Inner matrix M has dimensions (nFrames, jpj, jpi )
    Performs average passing through natural logarithm
    Mean  = exp(ln(values).mean() ) 
    At the moment works only for files on native mesh (253,733)
     
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

def getnextIndex(array,value):
    for i,val in enumerate(array):
        if val>value:
            break
    return i


    
    