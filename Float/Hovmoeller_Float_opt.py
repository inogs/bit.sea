from static.Float_opt_reader import Float_opt_reader
from commons.Timelist import TimeList, TimeInterval
from basins import V2 as OGS
import numpy as np
from instruments.statistics import mean_profile
from commons import genUserDateList as DL

def get_Hovmoeller(z, var, TL, region, Reader_Object):
    nLev = len(z)
    nFrames = TL.nTimes
    REQS=TL.getOwnList()
    HOV_MATRIX = np.zeros((nFrames,nLev), np.float32)*np.nan

    for iFrame in range(nFrames):
        req=REQS[iFrame]
        print req
        ProfileList = Reader_Object.Selector(var, req.time_interval, region)
        print len(ProfileList), "profiles"
        HOV_MATRIX[iFrame,:] = mean_profile(ProfileList, var, z)
    return HOV_MATRIX


z= np.arange(0,151,10)
N = Float_opt_reader()
var='chl'




TI = TimeInterval("20130501","20140501", '%Y%m%d')
monthly=DL.getTimeList(TI.start_time.strftime("%Y%m%d-%H:%M:%S"), TI.end_time.strftime("%Y%m%d-%H:%M:%S"), "months=1")
TL = TimeList(monthly)

#INPUTDIR="/marconi_work/OGS_dev_0/DEGRADATION_4_70/TEST_16/wrkdir/MODEL/AVE_FREQ_2/"
#TL = TimeList.fromfilenames(TI, INPUTDIR, "*nc", filtervar="N1p")



HOV_MATRIX = get_Hovmoeller(z, var, TL, OGS.ion, N)


import matplotlib.pyplot as pl
import numpy.ma as ma
Zm = ma.masked_where(np.isnan(HOV_MATRIX),HOV_MATRIX)

fig, ax = pl.subplots()
#xs,ys = np.meshgrid(TL.Timelist, z)
quadmesh = ax.pcolormesh(TL.Timelist, z, Zm.T, shading='flat')




#m = mean_profile(ProfileList,'chl',z)