import numpy as np
import scipy.io.netcdf as NC
import netcdf3




def TimeAverager3D(Filelist,weights,varname,mask):
    '''
    Performs a weighted average, working on a list of files.
    The weights are usually provided by TimeList.select() method.
    varname is a string
    mask is an common.mask.Mask object
    '''
    n=len(Filelist)
    jpk, jpj, jpi= mask.shape
    MSUM=np.zeros((jpk,jpj,jpi),np.float32)
    for t in range(n):
        filename=Filelist[t]
        #De      = DataExtractor(TheMask,filename,varname)
        #M = De.values
        M = netcdf3.read_3d_file(filename, varname)
        MSUM += M*weights[t]
    averaged = MSUM/weights.sum()
    return averaged

def TimeAverager2D(Filelist,weights,varname,mask):
    '''
    Performs a weighted average, working on a list of files.
    The weights are usually provided by TimeList.select() method.
    varname is a string
    mask is an common.mask.Mask object
    '''
    n=len(Filelist)
    _,jpj, jpi= mask.shape
    MSUM=np.zeros((jpj,jpi),np.float32)
    for t in range(n):
        filename=Filelist[t]
        #De      = DataExtractor(TheMask,filename,varname)
        #M = De.values
        ncIN = NC.netcdf_file(filename,"r")
        M = ncIN.variables[varname].data[:,:].copy()
        ncIN.close()
        MSUM += M*weights[t]
    averaged = MSUM/weights.sum()
    return averaged


def TimeAverager3D_2(Filelist,weights,varname,mask):
    '''
    Performs a weighted average of the square of a var (needed for STD on time dimension), working on a list of files.
    The weights are usually provided by TimeList.select() method.
    varname is a string
    mask is an common.mask.Mask object
    '''
    n=len(Filelist)
    jpk, jpj, jpi= mask.shape
    MSUM=np.zeros((jpk,jpj,jpi),np.float32)
    for t in range(n):
        filename=Filelist[t]
        #De      = DataExtractor(TheMask,filename,varname)
        #M = De.values
        M = netcdf3.read_3d_file(filename, varname)
        MSUM += M*M*weights[t]
    averaged = MSUM/weights.sum()
    return averaged
