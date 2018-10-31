from commons import netcdf4
from commons import timerequestors
from commons.Timelist import TimeInterval, TimeList
from Sat import SatManager
from postproc import masks
import numpy as np
import os,sys

maskSat = getattr(masks,"Mesh24")
jpi = maskSat.jpi
jpj = maskSat.jpj


MODELDIR="/pico/scratch/userexternal/gbolzon0/EOF-python/ORIG_INPUTS/CHL_SUP24/"
SATDIR  ="/pico/scratch/userexternal/gbolzon0/EOF-python/ORIG_INPUTS/CCI1km_Interp24_weekly/"
OUTPUTDIR="/pico/scratch/userexternal/gbolzon0/EOF-python/VarErr/"



TL_mod = TimeList.fromfilenames(None, MODELDIR, "*nc", prefix="chl.")
TL_sat = TimeList.fromfilenames(None, SATDIR, "*.nc", prefix="", dateformat="%Y%m%d")
suffix = os.path.basename(TL_sat.filelist[0])[8:]


for month in range(1,13):
    req = timerequestors.Clim_month(month)
    ii,w =TL_mod.select(req)
    
    sumDiff2 = np.zeros((jpj, jpi), np.float32)
    Npoints = np.zeros((jpj, jpi) , np.int32)

    model_filelist=[TL_mod.filelist[k] for k in ii if TL_mod.Timelist[k].month==month] # we consider strictly the month
    model_timelist=[TL_mod.Timelist[k] for k in ii if TL_mod.Timelist[k].month==month]
    
    for iFrame, modelfile in enumerate(model_filelist):
        modeltime=model_timelist[iFrame]
        CoupledList = TL_sat.couple_with([modeltime])
        sattime = CoupledList[0][0]
        assert modeltime == sattime + timerequestors.relativedelta(days=3.5)
        satfile = SATDIR + sattime.strftime("%Y%m%d") + suffix
        Mod = SatManager.readfromfile(modelfile,'lchlm')
        Sat = SatManager.readfromfile(satfile)
        diff = Sat - Mod
        cloudsLand = (np.isnan(Sat)) | (Sat > 1.e19) | (Sat<0)
        modelLand  = Mod == 1.e+20
        nodata     = cloudsLand | modelLand
        diff[nodata] = 0
        Npoints = Npoints + (~nodata).astype(np.int32)
        sumDiff2 += diff**2

    
    Npoints[Npoints==0] = 1 # just to divide
    varDiff = sumDiff2/Npoints
    nodata = varDiff==0
    varDiff[nodata] = -999    

    outfile="%svarErr.%02d.nc"  %( OUTPUTDIR, month)
    SatManager.dumpGenericNativefile(outfile, varDiff, "variance", maskSat)