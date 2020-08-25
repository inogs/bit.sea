from postproc import masks
import scipy.io.netcdf as NC
import numpy as np
import netCDF4

def mean(array):
    return array.mean()
def logmean(array):
    LOGmean= np.log(array).mean()
    return np.exp(LOGmean)

NativeMesh = masks.SatOrigMesh
OneKmMesh  = masks.SAT1km_mesh
V4 = masks.V4mesh
fillValue = -999.0

    
def readfromfile(filename,var='CHL'):
    '''
    returns CHL
    ''' 
    ncIN = netCDF4.Dataset(filename,'r')
    varObj = ncIN.variables[var]
    ndims = len(varObj.shape)
    if ndims==2: CHL_IN=np.array(varObj)
    if ndims==3: CHL_IN=np.array(varObj[0,:,:])
    del varObj
    ncIN.close()
    return CHL_IN

def readClimatology(filename):
    ncIN = netCDF4.Dataset(filename,'r')
    MEAN = np.array( ncIN.variables['Mean'])
    STD  = np.array( ncIN.variables['Std'])#(366, 253, 733)
    ncIN.close()
    return MEAN,STD

def readBathymetry(filename="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/MODIS/Bathy_MODIS.nc"):
    ncIN = netCDF4.Dataset(filename,'r')
    BATHY=np.array(ncIN.variables['bathy'])
    ncIN.close()
    return BATHY

def filterOnBathymetry(chl,Bathy,bathy_threshold=0):
    CHL=chl.copy()
    CHL[Bathy<bathy_threshold]=fillValue
    return CHL

def convertinV4format(CHL):
    '''
    Returns V4 formatted matrix having nans
    '''
    jpj,jpi = CHL.shape
    assert jpj == NativeMesh.jpj
    assert jpi == NativeMesh.jpi
    chl = np.ones((V4.jpj, V4.jpi) ,np.float32) * fillValue
    chl[1:,:] = CHL[0:-1:,11:]
    chl[chl==fillValue] = np.nan
    return chl

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
    BATHY=readBathymetry()
    chl16f = filterOnBathymetry(CHL_16, BATHY, bathy_threshold)
    # this is just an hint of how to develop/simplify the next part
    
    PREVLON = np.array( [getnextIndex(NativeMesh.lon, V1.lon[i]) for i in range(V1.jpi)], np.int32) -1
    PREVLAT = np.array( [getnextIndex(NativeMesh.lat, V1.lat[i]) for i in range(V1.jpj)], np.int32) -1

    CHL_8 = np.ones((V1.jpj,V1.jpi),np.float32) * fillValue

    for i in range(V1.jpi):
        ji = PREVLON[i]
        for j in range(V1.jpj):
            jj = PREVLAT[j] 

            riq =CHL_16[jj:jj+2,ji:ji+2]
            bat = BATHY[jj:jj+2,ji:ji+2]

            batgoods = bat > bathy_threshold
            condition = False

            if batgoods.sum() > 3:
                condition = True
            if batgoods.sum() == 2:
                diag1 = batgoods.diagonal()
                diag2 = np.array([batgoods[0,1], batgoods[1,0]],np.bool)
                if (diag1.all()) | (diag2.all()) :
                    condition=True

            if condition:
                selected = riq[batgoods]
                writeable = selected>fillValue
                if np.any(writeable):
                    CHL_8[j,i] = selected[writeable].mean()
    return CHL_8

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

def dumpGenericNativefile(filename, M, varname='KD490', mesh=masks.KD490mesh):
    '''
    Used in sat check
    '''
    
    ncOUT  = NC.netcdf_file(filename,'w')
    ncOUT.createDimension('time', 1)
    ncOUT.createDimension('depth',1)
    ncOUT.createDimension('lon',mesh.jpi)
    ncOUT.createDimension('lat',mesh.jpj)
    
    ncvar = ncOUT.createVariable(varname, 'f', ('time','lat','lon'))
    chl = np.zeros((1,mesh.jpj,mesh.jpi),np.float32)
    chl[0,:,:] = M
    ncvar[:] = chl
    setattr(ncvar, 'missing_value', fillValue)
    setattr(ncvar, '_FillValue', fillValue)
    ncvar = ncOUT.createVariable('depth','f',('depth',))
    ncvar[:] = 1.47210180759
    ncvar = ncOUT.createVariable('time','f',('time',))
    ncvar[:] = 1.0
    ncvar = ncOUT.createVariable('lat','f',('lat',))
    ncvar[:] = mesh.lat
    ncvar = ncOUT.createVariable('lon','f',('lon',))
    ncvar[:] = mesh.lon
    ncOUT.close()



def dumpV4file(filename,CHL,varname='lchlm'):
    '''
    Second argument CHL is a matrix in the native format
    '''
    ncOUT  = NC.netcdf_file(filename,'w')
    ncOUT.createDimension('time', 1)
    ncOUT.createDimension('depth',1)
    ncOUT.createDimension('lon',V4.jpi)
    ncOUT.createDimension('lat',V4.jpj)
    
    ncvar = ncOUT.createVariable(varname, 'f', ('lat','lon'))
    CHLv4 = convertinV4format(CHL)
    CHLv4[np.isnan(CHLv4)] = fillValue
    ncvar[:] = CHLv4
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


def dump_simple_V4file(filename,CHL,varname='lchlm'):
    '''
    Second argument CHL is a matrix in the V4 format
    '''
    ncOUT  = NC.netcdf_file(filename,'w')
    ncOUT.createDimension('time', 1)
    ncOUT.createDimension('depth',1)
    ncOUT.createDimension('lon',V4.jpi)
    ncOUT.createDimension('lat',V4.jpj)
    
    ncvar = ncOUT.createVariable(varname, 'f', ('lat','lon'))
    ncvar[:] = CHL
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
    

    '''
    _,jpj,jpi = M.shape
    CHL_OUT = np.ones((jpj,jpi),np.float32) * fillValue

    maskM = M==fillValue
    Mnan = M.copy()
    Mnan[maskM] = np.nan
    LOGm = np.nanmean(np.log(Mnan),0)
    CHL_OUT = np.exp(LOGm)
    CHL_OUT[np.isnan(CHL_OUT)] = fillValue
    return CHL_OUT

def averager(M):
    '''
    Inner matrix M has dimensions (nFrames, jpj, jpi )
    Performs average on present values (avoiding fillvalues)

    '''
    _,jpj,jpi = M.shape
    #assert jpj == NativeMesh.jpj
    #assert jpi == NativeMesh.jpi
    CHL_OUT = np.ones((jpj,jpi),np.float32) * fillValue

#    for i in range(jpi):
#        for j in range(jpj):
#            l = M[:,j,i]
#            goodValues = l != fillValue
#            if np.any(goodValues):
#                count = goodValues.sum()
#                CHL_OUT[j,i] = l[goodValues].sum() / count
    maskM = M==fillValue
    Mnan = M.copy()
    Mnan[maskM] = np.nan
    CHL_OUT = np.nanmean(Mnan,0)
    CHL_OUT[np.isnan(CHL_OUT)] = fillValue
    return CHL_OUT


def WeightedLogAverager(M, w):
    '''
    Inner matrix M has dimensions (nFrames, jpj, jpi ), while w has dimension (nFrames)
    Performs a log weighted average on present values (avoiding fillvalues)

    '''
    _,jpj,jpi = M.shape
    #assert jpj == NativeMesh.jpj
    #assert jpi == NativeMesh.jpi
    CHL_OUT = np.ones((jpj,jpi),np.float32) * fillValue

    nFrames = M.shape[0]
    maskM = (M==fillValue)==False
    maskNotallNan = np.sum(maskM,0)>0
    nP = np.sum(maskNotallNan)
    Mmasked = np.zeros((nFrames,nP))
    wwM = np.zeros_like(Mmasked)
    for iframe in range(nFrames):
        wwM[iframe,:] = w[iframe]
        Mmasked[iframe,:] = M[iframe,maskNotallNan]
    noValidM = Mmasked==fillValue
    Mmasked[noValidM] = np.nan
    wwM[noValidM] = np.nan
    countM = np.nansum(wwM,0)
    chlm = np.nansum(wwM*np.log(Mmasked),0) / countM
    CHL_OUT[maskNotallNan] = np.exp(chlm)

    return CHL_OUT


def WeightedAverager(M, w):
    '''
    Inner matrix M has dimensions (nFrames, jpj, jpi ), while w has dimension (nFrames)
    Performs a weighted average on present values (avoiding fillvalues)

    '''
    _,jpj,jpi = M.shape
    #assert jpj == NativeMesh.jpj
    #assert jpi == NativeMesh.jpi
    CHL_OUT = np.ones((jpj,jpi),np.float32) * fillValue

    nFrames = M.shape[0]
    maskM = (M==fillValue)==False
    maskNotallNan = np.sum(maskM,0)>0
    nP = np.sum(maskNotallNan)
    Mmasked = np.zeros((nFrames,nP))
    wwM = np.zeros_like(Mmasked)
    for iframe in range(nFrames):
        wwM[iframe,:] = w[iframe]
        Mmasked[iframe,:] = M[iframe,maskNotallNan]
    noValidM = Mmasked==fillValue
    Mmasked[noValidM] = np.nan
    wwM[noValidM] = np.nan
    countM = np.nansum(wwM,0)
    chlm = np.nansum(wwM*Mmasked,0) / countM
    CHL_OUT[maskNotallNan] = chlm

    return CHL_OUT



def getnextIndex(array,value):
    for i,val in enumerate(array):
        if val>value:
            break
    return i


    
    
