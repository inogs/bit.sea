from static.Float_opt_reader import Float_opt_reader
from commons.Timelist import TimeList, TimeInterval
from basins import V2 as OGS
import numpy as np


z= np.arange(0,150,10)
N = Float_opt_reader()
nLev = len(z)
var='chl'


def mean_profile(ProfileList, var, depth):
    '''
    Calculates the mean profile of Profilelist
    Arguments:
    * ProfileList * a list of Profile Objects
    * var         * string, the var name
    * depth       * numpy array, the depth of the output

    Returns:
    * MEAN * numpy array corresponding to depth where profiles of measurements are interpolated
    '''
    nProfiles = len(ProfileList) 
    nLev = len(z)
    
    M = np.zeros((nProfiles, nLev), np.float32)*np.nan
    MEAN = np.zeros((nLev,),np.float32) * np.nan
    
    for ip, p in enumerate(ProfileList):
        Pres, Profile, Qc = p.read(var)
        M[ip,:] = np.interp(depth, Pres, Profile, right=np.nan)
    for iLev in range(nLev):
        level= M[:,iLev]
        good = ~np.isnan(level)
        if good.sum() > 2:
            MEAN[iLev] = level[good].mean()
    return MEAN


if __name__=="__main__":
    TI = TimeInterval("20130501","20140501", '%Y%m%d')
    INPUTDIR="/marconi_work/OGS_dev_0/DEGRADATION_4_70/TEST_16/wrkdir/MODEL/AVE_FREQ_2/"
    TL = TimeList.fromfilenames(TI, INPUTDIR, "*nc", filtervar="N1p")
    nFrames = TL.nTimes
    REQS=TL.getOwnList()
    HOV_MATRIX = np.zeros((nFrames,nLev), np.float32)*np.nan

    for iFrame in range(nFrames):
        req=REQS[iFrame]
        print req
        ProfileList = N.Selector(var, req.time_interval, OGS.ion)
        print len(ProfileList), "profiles"
        HOV_MATRIX[iFrame,:] = mean_profile(ProfileList, var, z)



    import pylab as pl
    import numpy.ma as ma
    Zm = ma.masked_where(np.isnan(HOV_MATRIX),HOV_MATRIX)

    fig, ax = pl.subplots()
    #xs,ys = np.meshgrid(TL.Timelist, z)
    quadmesh = ax.pcolormesh(TL.Timelist, z, Zm.T, shading='flat')




#m = mean_profile(ProfileList,'chl',z)