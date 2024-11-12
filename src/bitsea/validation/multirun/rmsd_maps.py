import matplotlib
matplotlib.use('Qt5Agg')
import numpy as np
from bitsea.commons.mask import Mask
from bitsea.commons.dataextractor import DataExtractor
from bitsea.commons.Timelist import TimeList, TimeInterval
import bitsea.matchup.matchup as matchup
from bitsea.commons import netcdf4
import pylab as pl

TI = TimeInterval("20190201","20190501"); season='winter'
#TI = TimeInterval("20190701","20191001"); season='summer'


TheMask = Mask.from_file('/g100_work/OGS_devC/gbolzon/BI-HOURLY/SETUP/meshmask.nc')
mask2d=TheMask.mask_at_level(0)
jpk,jpj,jpi = TheMask.shape
MODEL_24="/g100_scratch/userexternal/gbolzon0/BI-HOURLY/24H_orig/wrkdir/POSTPROC/output/AVE_FREQ_1/CHL_SUP/"
MODEL_2H="/g100_scratch/userexternal/gbolzon0/BI-HOURLY/2H/wrkdir/POSTPROC/output/AVE_FREQ_1/CHL_SUP/"
outfile_24=season + "_RMSE_24H_orig.nc"
outfile_2H=season + "RMSE_2H.nc"
REF_DIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE_V9C/SAT/CHL/DT/DAILY/CHECKED_24/"
dateformat ="%Y%m%d"



def get_timeseries(filelist,varname):
    n = len(filelist)
    M = np.zeros((n,jpj,jpi), np.float32)
    for i, filename in enumerate(filelist):
        M[i,:] = DataExtractor(TheMask,filename,varname,dimvar=2).values
    return M


def get_rmsmap(Model,Ref):
    '''
    Model and ref are 3D ndrarrays [nFrames, jpj,jpi]
    '''
    RMSmap=np.ones((jpj,jpi),np.float32)*1.e+20
    NUMB = np.ones((jpj,jpi), np.float32)*1.e+20

    for ji in range(jpi):
        for jj in range(jpj):
            if not mask2d[jj,ji] : continue

            good = Ref[:,jj,ji] > -999
            n = good.sum()
            if n > 0:
                model = Model[good,jj,ji]
                obs   = Ref[good,jj,ji]
                M = matchup.matchup(model, obs)
                RMSmap[jj,ji] = M.RMSE()
                NUMB[jj,ji] = n
    return RMSmap, NUMB


sat_TL   = TimeList.fromfilenames(TI, REF_DIR  ,"*.nc", prefix="", dateformat=dateformat)
Mod2H_TL = TimeList.fromfilenames(TI, MODEL_2H,"*.nc", filtervar="P_l")
Mod24HTL = TimeList.fromfilenames(TI, MODEL_24,"*.nc", filtervar="P_l")

SAT =   get_timeseries(sat_TL.filelist, 'CHL')
MOD2H = get_timeseries(Mod2H_TL.filelist, 'lchlm')
MOD24 = get_timeseries(Mod24HTL.filelist, 'lchlm')


RMS_24, NUMB = get_rmsmap(MOD24, SAT)
RMS_2H, NUMB = get_rmsmap(MOD2H, SAT)

netcdf4.write_2d_file(RMS_24, 'RMS', outfile_24, TheMask)
netcdf4.write_2d_file(NUMB, 'NUMB', outfile_24, TheMask)
netcdf4.write_2d_file(RMS_2H, 'RMS', outfile_2H, TheMask)
netcdf4.write_2d_file(NUMB, 'NUMB', outfile_2H, TheMask)




def plotta_fascio(ji,jj):
    '''
    Plots timeseries of
    '''
    Largh_fascio=10 #points
    delta = int(Largh_fascio/2)
    fig,ax=pl.subplots()
    for i in range(ji-delta,ji+delta):
        for j in range(jj-delta,jj+delta):
            ax.plot(MOD24[:,j,i],'b')
            ax.plot(MOD2H[:,j,i],'r')
            good = SAT[:,j,i]>0
            t = np.nonzero(good)[0]
            ax.plot(t,SAT[t,j,i],'c.')

    # just to add legend
    ax.plot(MOD24[:,jj,ji],'b',label='24h')
    ax.plot(MOD2H[:,jj,ji],'r',label='2h')
    good = SAT[:,jj,ji]>0
    t = np.nonzero(good)[0]
    ax.plot(t, SAT[t,jj,ji],'c.',label='sat')
    ax.legend()
    ax.grid()
    fig.show()



#Caso invernale, nwm
ji,jj=365,279

plotta_fascio(365, 279)

ji,jj=635,275  # Sud-Adriatico

ji,jj = 514,346 #Po
ji,jj = 987,35# Nile
ji,jj = 341,196 # swm







