import numpy as np
from commons.mask import Mask
from commons.dataextractor import DataExtractor
from commons.Timelist import TimeList, TimeInterval
import matchup.matchup as matchup
from commons import netcdf4
TI = TimeInterval("20190201","20190501")


TheMask = Mask('filename')
mask2d=TheMask.mask_at_level(0)
jpk,jpj,jpi = TheMask.shape
MODEL_DIR=""
REF_DIR=""
dateformat ="%Y%m%d"
sat_TL   = TimeList.fromfilenames(TI, REF_DIR  ,"*.nc", prefix="", dateformat=dateformat)
model_TL = TimeList.fromfilenames(TI, MODEL_DIR,"*.nc", filtervar="P_l")

nFrames=model_TL.nTimes
Model = np.zeros((nFrames,jpj,jpi), np.float32)
Ref   = np.zeros((nFrames,jpj,jpi), np.float32)


for i, filename in enumerate(model_TL):
    D = DataExtractor(TheMask,filename,'P_l').values[0,:]
    Model[i,:] =  D

for i, filename in enumerate(sat_TL):
    D = DataExtractor(TheMask,filename,'P_l').values[0,:]
    Ref[i,:] =  D



RMSmap=np.ones((jpj,jpi),np.float32)*1.e+20
NUMB = np.zeros((jpj,jpi), np.int32)

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

netcdf4.write_2d_file(RMSmap, 'RMS', 'out.nc', TheMask)
