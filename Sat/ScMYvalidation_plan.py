# OUTPUTS
# BGC_CLASS4_CHL_RMS_SURF_BASIN
# of CMEMS-Med-biogeochemistry-S.0.pdf

from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
import commons.IOnames as IOnames
import numpy as np
import SatManager as Sat
import matchup.matchup as matchup
from commons.dataextractor import DataExtractor
from layer_integral.mapbuilder import MapBuilder
from commons.mask import Mask
from commons.submask import SubMask
from basins import OGS
from commons.layer import Layer

def weighted_mean(Conc, Weight):
    
    Weight_sum      = Weight.sum()
    Mass            = (Conc * Weight).sum()
    Weighted_Mean   = Mass/Weight_sum
    return Weighted_Mean 

TheMask=Mask('/pico/home/usera07ogs/a07ogs00/OPA/V4/etc/static-data/MED1672_cut/MASK/meshmask.nc')
MODEL_DIR="/pico/scratch/userexternal/gbolzon0/RA_CARBO/RA/wrkdir/POSTPROC/output/AVE_FREQ_2/TMP/"
REF_DIR  = "/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/CCI/MONTHLY/"

Timestart="19990101"
Time__end="20120731"

TI    = TimeInterval(Timestart,Time__end,"%Y%m%d")
IonamesFile = '../postproc/IOnames_sat.xml'
IOname = IOnames.IOnames(IonamesFile)
sat_TL = TimeList.fromfilenames(TI, REF_DIR,"*.nc",IonamesFile)

IonamesFile = '../postproc/IOnames.xml'
model_TL = TimeList.fromfilenames(TI, MODEL_DIR,"*.nc",IonamesFile)


nFrames = model_TL.nTimes
nSUB = len(OGS.P.basin_list)

jpk,jpj,jpi =TheMask.shape
dtype = [(sub.name, np.bool) for sub in OGS.P]
SUB = np.zeros((jpj,jpi),dtype=dtype)
for sub in OGS.P:
    sbmask         = SubMask(sub,maskobject=TheMask).mask
    SUB[sub.name]  = sbmask[0,:,:]

mask200_2D = TheMask.mask_at_level(200.0)

BGC_CLASS4_CHL_RMS_SURF_BASIN      = np.zeros((nFrames,nSUB),np.float32)
BGC_CLASS4_CHL_BIAS_SURF_BASIN     = np.zeros((nFrames,nSUB),np.float32)
MODEL_MEAN                         = np.zeros((nFrames,nSUB),np.float32)
SAT___MEAN                         = np.zeros((nFrames,nSUB),np.float32)
# This is the surface layer choosen to match satellite chl data
surf_layer = Layer(0,10)

for itime, modeltime in enumerate(model_TL.Timelist):
    print modeltime
    CoupledList = sat_TL.couple_with([modeltime])
    sattime = CoupledList[0][0]
    satfile = REF_DIR + sattime.strftime(IOname.Input.dateformat) + IOname.Output.suffix + ".nc"
    modfile = model_TL.filelist[itime]
     
    De         = DataExtractor(TheMask,filename=modfile, varname='P_i')
    Model      = MapBuilder.get_layer_average(De, surf_layer)
    #ncIN = NC.netcdf_file(modfile,'r')
    #Model = ncIN.variables['P_i'].data[0,0,:,:].copy()#.astype(np.float64)
    #Model = ncIN.variables['lchlm'].data.copy()
    #ncIN.close()
    
    Sat16 = Sat.readfromfile(satfile,var='lchlm') #.astype(np.float64)
    
    
    cloudsLand = (np.isnan(Sat16)) | (Sat16 > 1.e19) | (Sat16<0)
    modelLand  = np.isnan(Model) #lands are nan
    nodata     = cloudsLand | modelLand
    selection = ~nodata & mask200_2D 
    M = matchup.matchup(Model[selection], Sat16[selection])
    
    for isub, sub in enumerate(OGS.P):
        selection = SUB[sub.name] & (~nodata) & mask200_2D
        M = matchup.matchup(Model[selection], Sat16[selection])
        BGC_CLASS4_CHL_RMS_SURF_BASIN[itime,isub]  = M.RMSE()
        BGC_CLASS4_CHL_BIAS_SURF_BASIN[itime,isub] = M.bias()
        weight = TheMask.area[selection]
        MODEL_MEAN[itime,isub] = weighted_mean( M.Model,weight)
        SAT___MEAN[itime,isub] = weighted_mean( M.Ref,  weight)
    
BGC_CLASS4_CHL_EAN_RMS_SURF_BASIN  = BGC_CLASS4_CHL_RMS_SURF_BASIN.mean(axis=0)
BGC_CLASS4_CHL_EAN_BIAS_SURF_BASIN = BGC_CLASS4_CHL_BIAS_SURF_BASIN.mean(axis=0)


import pickle
LIST=[model_TL.Timelist, BGC_CLASS4_CHL_RMS_SURF_BASIN, BGC_CLASS4_CHL_BIAS_SURF_BASIN,MODEL_MEAN, SAT___MEAN]
fid = open('export_data_ScMYValidation_plan.pkl','wb')
pickle.dump(LIST, fid)
fid.close()