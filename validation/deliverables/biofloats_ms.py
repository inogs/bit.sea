import argparse

from basins import V2 as OGS
from instruments import lovbio_float as bio_float
from instruments.var_conversions import LOVFLOATVARS
from instruments.matchup_manager import Matchup_Manager
from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
from commons.mask import Mask
from commons.layer import Layer
import numpy as np
from matchup.statistics import matchup
import datetime
import scipy.io.netcdf as NC
from commons.utils import addsep
from basins.region import Rectangle
from profiler import ALL_PROFILES, TL, BASEDIR

maskfile="/pico/scratch/userexternal/gbolzon0/eas_v12/eas_v12_8/wrkdir/MODEL/meshmask.nc"
TheMask  = Mask(maskfile)



LAYERLIST=[Layer(0,10), Layer(10,30), Layer(30,60), Layer(60,100), Layer(100,150), Layer(150,300), Layer(300,600), Layer(600,1000)]
VARLIST = ['P_l','N3n','O2o']
read_adjusted = [True,False,False]
nSub   = len(OGS.NRT3.basin_list)
nDepth = len(LAYERLIST)
nVar   = len(VARLIST)


nFrames = TL.nTimes
BIAS    = np.zeros((nVar,nFrames,nSub,nDepth), np.float32)
RMSE    = np.zeros((nVar,nFrames,nSub,nDepth), np.float32)
NPOINTS = np.zeros((nVar,nFrames, nSub,nDepth), np.int32)


M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)



for iFrame, req in enumerate(TL.getOwnList()):
    print req.time_interval
    for ivar, var in enumerate(VARLIST):
        print var
        for isub, sub in enumerate(OGS.NRT3):
            Profilelist = bio_float.FloatSelector(LOVFLOATVARS[var], req.time_interval, sub)
            nProfiles = len(Profilelist)
            print sub.name, nProfiles
            Matchup_object_list=[]
            for ip in range(nProfiles):
                floatmatchup =  M.getMatchups([Profilelist[ip]], TheMask.zlevels, var, read_adjusted=read_adjusted[ivar])
                Matchup_object_list.append(floatmatchup)
    
            for ilayer, layer in enumerate(LAYERLIST):
                MODEL_LAYER_MEAN = [] # one value for each suitable profile in (subbasin, layer)
                REF_LAYER_MEAN   = []
                for floatmatchup in Matchup_object_list:
                    m_layer = floatmatchup.subset(layer)
                    #print ilayer, m_layer.number()
                    if m_layer.number() > 0:
                        REF_LAYER_MEAN.append(m_layer.Ref.mean())
                        MODEL_LAYER_MEAN.append(m_layer.Model.mean())
    
                NPOINTS[ivar, iFrame, isub, ilayer] = len(MODEL_LAYER_MEAN)
                if len(MODEL_LAYER_MEAN) > 0:
                    M_LAYER = matchup(np.array(MODEL_LAYER_MEAN), np.array(REF_LAYER_MEAN))
                    BIAS[ivar, iFrame, isub, ilayer] = M_LAYER.bias()
                    RMSE[ivar, iFrame, isub, ilayer] = M_LAYER.RMSE()




outfile = "pippo.nc"
ncOUT = NC.netcdf_file(outfile,'w') #manca l'array times
ncOUT.createDimension('time', nFrames)
ncOUT.createDimension('var', nVar)
ncOUT.createDimension('sub', nSub)
ncOUT.createDimension('depth',nDepth)
s=''
for var in VARLIST: s= s+var + ","
setattr(ncOUT, 'varlist',s[:-1])
s='';
for sub in OGS.NRT3: s =s+sub.name + ","
setattr(ncOUT,'sublist',s[:-1])
s='';
for layer in LAYERLIST: s =s+layer.string() + ","
setattr(ncOUT,'layerlist',s[:-1])

ncvar=ncOUT.createVariable('bias', 'f', ('var','time', 'sub','depth'))
ncvar[:] = BIAS
ncvar=ncOUT.createVariable('rmse', 'f', ('var','time','sub','depth'))
ncvar[:] = RMSE
ncvar=ncOUT.createVariable('npoints', 'i', ('var','time','sub','depth'))
ncvar[:] = NPOINTS

ncOUT.close()
