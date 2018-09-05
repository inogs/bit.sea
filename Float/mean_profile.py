from static.Float_opt_reader import Float_opt_reader
from commons.Timelist import TimeList, TimeInterval
from basins import V2 as OGS
import numpy as np

DATESTART = '20150101-00:00:00'
DATE__END = '20150201-00:00:00'

TI = TimeInterval(DATESTART,DATE__END, '%Y%m%d-%H:%M:%S')
z= np.arange(0,200,10)
N = Float_opt_reader()
nLev = len(z)
var='chl'
ProfileList = N.Selector('chl',TI,OGS.ion)


def mean_profile(ProfileList, var, depth):
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



TI = TimeInterval("20130101","20160101", '%Y%m%d')
INPUTDIR="/marconi_work/OGS_dev_0/DEGRADATION_4_70/TEST_16/wrkdir/MODEL/AVE_FREQ_2/"
TL = TimeList.fromfilenames(TI, INPUTDIR, "*nc", filtervar="N1p")
nFrames = TL.nTimes
REQS=TL.getOwnList()
HOV_MATRIX = np.zeros((nFrames,nLev), np.float32)*np.nan

for iFrame in range(nFrames):
    req=REQS[iFrame]
    print req
    ProfileList = N.Selector(var, req.time_interval, OGS.ion)
    HOV_MATRIX[iFrame,:] = mean_profile(ProfileList, var, z)

import pylab as pl
fig, ax = pl.subplots()
quadmesh = ax.pcolormesh(TL.Timelist, z, HOV_MATRIX,shading='gouraud')




#m = mean_profile(ProfileList,'chl',z)