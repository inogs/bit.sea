import numpy as np
from commons import netcdf3
from commons.dataextractor import DataExtractor




def TimeAverager3D(Filelist,weights,varname,mask,accept_neg=True):
    '''
    Performs a weighted average, working on a list of files.
    The weights are usually provided by TimeList.select() method.
    varname is a string
    mask is an common.mask.Mask object
    if accept_neg is False, it modifies negative values in zeros
    '''
    n=len(Filelist)
    jpk, jpj, jpi= mask.shape
    MSUM=np.zeros((jpk,jpj,jpi),np.float32)
    for t in range(n):
        filename=Filelist[t]
        De      = DataExtractor(mask,filename,varname)
        M = De.values
        #M = netcdf3.read_3d_file(filename, varname)
        if not accept_neg:
             M[M<0]=0
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
        De      = DataExtractor(mask,filename,varname,dimvar=2)
        M = De.values

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


def TimeAverager3D_std(Filelist,weights,varname,mask):
    '''
    Performs a weighted average and std dev, working on a list of files.
    This function uses Welford's one loop formula for the 
    computation of std_dev computation.
    The weights are usually provided by TimeList.select() method.
    varname is a string
    mask is an common.mask.Mask object
    '''
    n             = len(Filelist)
    jpk, jpj, jpi = mask.shape
    average       = np.zeros((jpk,jpj,jpi),np.float32)
    std_dev       = np.zeros((jpk,jpj,jpi),np.float32)
    sum_weight    = 0.
    for t in range(n):
        # Extract the data
        filename = Filelist[t]
        De       = DataExtractor(mask,filename,varname)
        M        = De.values
        # Sum the weight
        sum_weight += weights[t]
        # Compute mean and std using Welford algorithm
        average_old = average.copy()
        average    += weights[t]/sum_weight*(M-average)
        std_dev    += weights[t]*(M-average_old)*(M-average)
    std_dev = np.sqrt(std_dev/sum_weight)
    return average, std_dev


def TimeAverager2D_std(Filelist,weights,varname,mask):
    '''
    Performs a weighted average and std dev, working on a list of files.
    This function uses Welford's one loop formula for the 
    computation of std_dev computation.
    The weights are usually provided by TimeList.select() method.
    varname is a string
    mask is an common.mask.Mask object
    '''
    n             = len(Filelist)
    _,jpj, jpi    = mask.shape
    average       = np.zeros((jpj,jpi),np.float32)
    std_dev       = np.zeros((jpj,jpi),np.float32)
    sum_weight    = 0.
    for t in range(n):
        # Extract the data
        filename = Filelist[t]
        De       = DataExtractor(mask,filename,varname,dimvar=2)
        M        = De.values
        # Sum the weight
        sum_weight += weights[t]
        # Compute mean and std using Welford algorithm
        average_old = average.copy()
        average    += weights[t]/sum_weight*(M-average)
        std_dev    += weights[t]*(M-average_old)*(M-average)
    std_dev = np.sqrt(std_dev/sum_weight)
    return average, std_dev
