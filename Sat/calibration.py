# OUTPUTS
# BGC_CLASS4_CHL_EAN_RMS_SURF_BASIN
# BGC_CLASS4_CHL_EAN_BIAS_SURF_BASIN
# BGC_CLASS4_CHL_RMS_SURF_BASIN
# BGC_CLASS4_CHL_BIAS_SURF_BASIN
# of CMEMS-Med-biogeochemistry-ScCP-1.0.pdf

from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
import commons.IOnames as IOnames
import numpy as np
import SatManager as Sat
import matchup.matchup as matchup
from postproc.maskload import *
import scipy.io.netcdf as NC



MODEL_DIR="/pico/scratch/userexternal/gbolzon0/Carbonatic-17/wrkdir/POSTPROC/output/AVE_FREQ_1/CHL_SUP/"
REF_DIR  = "/pico/scratch/userexternal/gbolzon0/Carbonatic-01/SAT16/"

Timestart="20140404"
Time__end="20150701"

TI    = TimeInterval(Timestart,Time__end,"%Y%m%d")
IonamesFile = '../postproc/IOnames_sat.xml'
IOname = IOnames.IOnames(IonamesFile)
sat_TL = TimeList.fromfilenames(TI, REF_DIR,"*.nc",IonamesFile)

IonamesFile = '../postproc/IOnames.xml'
model_TL = TimeList.fromfilenames(TI, MODEL_DIR,"*.nc",IonamesFile)


nFrames = model_TL.nTimes
nSUB = len(SUBlist)

BGC_CLASS4_CHL_RMS_SURF_BASIN      = np.zeros((nFrames,nSUB),np.float32)
BGC_CLASS4_CHL_BIAS_SURF_BASIN     = np.zeros((nFrames,nSUB),np.float32)


for itime, modeltime in enumerate(model_TL.Timelist):
    print modeltime
    CoupledList = sat_TL.couple_with([modeltime])
    sattime = CoupledList[0][0]
    satfile = REF_DIR + sattime.strftime(IOname.Input.dateformat) + IOname.Output.suffix + ".nc"
    modfile = model_TL.filelist[itime]
     
    ncIN = NC.netcdf_file(modfile,'r')
    #Model = ncIN.variables['P_i'].data[0,0,:,:].copy()#.astype(np.float64)
    Model = ncIN.variables['lchlm'].data.copy()
    ncIN.close()
    
    Sat16 = Sat.readfromfile(satfile,var='lchlm') #.astype(np.float64)
    
    
    cloudsLand = (np.isnan(Sat16)) | (Sat16 > 1.e19)
    modelLand  = Model > 1.0e+19
    nodata     = cloudsLand | modelLand
    selection = ~nodata & mask200_2D 
    M = matchup.matchup(Model[selection], Sat16[selection])
    
    for isub, sub in enumerate(SUBlist):
        selection = SUB[sub][0,:,:] & (~nodata) & mask200_2D
        M = matchup.matchup(Model[selection], Sat16[selection])
        BGC_CLASS4_CHL_RMS_SURF_BASIN[itime,isub]  = M.RMSE()
        BGC_CLASS4_CHL_BIAS_SURF_BASIN[itime,isub] = M.bias()
    
BGC_CLASS4_CHL_EAN_RMS_SURF_BASIN  = BGC_CLASS4_CHL_RMS_SURF_BASIN.mean(axis=0)
BGC_CLASS4_CHL_EAN_BIAS_SURF_BASIN = BGC_CLASS4_CHL_BIAS_SURF_BASIN.mean(axis=0)