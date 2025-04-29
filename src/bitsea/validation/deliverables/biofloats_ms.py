import argparse
from bitsea.utilities.argparse_types import existing_dir_path, existing_file_path, generic_path
def argument():
    parser = argparse.ArgumentParser(description = '''
    Produces a single file, containing bias, rmse and number of measurements for each subbasin and layer
    for chlorophyll, nitrate and oxygen.
    In this approach we define the measurement as mean on layer of the float profile values.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--maskfile', '-m',
                                type = existing_file_path,
                                required = False)
    parser.add_argument(   '--basedir', '-b',
                                type = existing_dir_path,
                                required = True,
                                help = """ PROFILATORE dir, already generated""")
    parser.add_argument(   '--outfile', '-o',
                                type = generic_path,
                                required = True,
                                help = "")

    return parser.parse_args()

args = argument()

from bitsea.basins import V2 as OGS
from bitsea.instruments import superfloat as bio_float
from bitsea.instruments.var_conversions import FLOATVARS
from bitsea.instruments.matchup_manager import Matchup_Manager
from bitsea.commons.Timelist import TimeList
from bitsea.commons.time_interval import TimeInterval
from bitsea.commons.mask import Mask
from bitsea.commons.layer import Layer
import numpy as np
from bitsea.matchup.statistics import matchup
import datetime
import netCDF4
from pathlib import Path

from bitsea.basins.region import Rectangle
from bitsea.instruments import check
BASEDIR = args.basedir
TL = TimeList.fromfilenames(None, BASEDIR / "PROFILES/","ave*.nc")
deltaT= datetime.timedelta(hours=12)
TI = TimeInterval.fromdatetimes(TL.Timelist[0] - deltaT, TL.Timelist[-1] + deltaT)
ALL_PROFILES = bio_float.FloatSelector(None, TI, Rectangle(-6,36,30,46))
M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)

Check_obj_nitrate = check.check(Path(""), verboselevel=0)
Check_obj_chl     = check.check(Path(""), verboselevel=0)
Check_obj_PhytoC  = check.check(Path(""), verboselevel=0)

TheMask = Mask.from_file(args.maskfile)




LAYERLIST=[Layer(0,10), Layer(10,30), Layer(30,60), Layer(60,100), Layer(100,150), Layer(150,300), Layer(300,600), Layer(600,1000)]
VARLIST = ['P_l','N3n','O2o','P_c']

#read_adjusted = [True,True,True,True]
extrap = [True,False,False,False,False]
nSub   = len(OGS.NRT3.basin_list)
nDepth = len(LAYERLIST)
nVar   = len(VARLIST)


WEEKLY=TL.getWeeklyList(5)
nFrames = len(WEEKLY)
BIAS    = np.zeros((nVar,nFrames,nSub,nDepth), np.float32)*np.nan
RMSE    = np.zeros((nVar,nFrames,nSub,nDepth), np.float32)*np.nan
NPOINTS = np.zeros((nVar,nFrames, nSub,nDepth), np.int32)*np.nan

MOD = np.zeros((nVar,nFrames,nSub,nDepth), np.float32)*np.nan
REF = np.zeros((nVar,nFrames,nSub,nDepth), np.float32)*np.nan

M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)


#LISTcheck = []
#LISTall = []

for iFrame, req in enumerate(WEEKLY):
    if req.time_interval.start_time < TI.start_time : req.time_interval.start_time = TI.start_time
    if req.time_interval.end_time   > TI.end_time   : req.time_interval.end_time   = TI.end_time
    print (req, flush=True)
    for ivar, var in enumerate(VARLIST):
        if var == "N3n": Check_obj = Check_obj_nitrate
        if var == "P_l": Check_obj = Check_obj_chl
        if var == "O2o": Check_obj = None
        if var == "P_c": Check_obj = Check_obj_PhytoC
        if var == "POC": Check_obj = None

        for isub, sub in enumerate(OGS.NRT3):
            Profilelist_raw = bio_float.FloatSelector(FLOATVARS[var], req.time_interval, sub)
            Profilelist = bio_float.remove_bad_sensors(Profilelist_raw,FLOATVARS[var])
            nProfiles = len(Profilelist)
            Matchup_object_list=[]
            for ip in range(nProfiles):
                floatmatchup =  M.getMatchups2([Profilelist[ip]], TheMask.zlevels, var, interpolation_on_Float=False,checkobj=Check_obj, extrapolation=extrap[ivar])
#                LISTall.append(req.string + var + sub.name)
#                if all(floatmatchup.CheckReports)==None:
#                    LISTcheck.append(req.string + var + sub.name)

                Matchup_object_list.append(floatmatchup)
    
            for ilayer, layer in enumerate(LAYERLIST):
                MODEL_LAYER_MEAN = [] # one value for each suitable profile in (subbasin, layer)
                REF_LAYER_MEAN   = []
                for floatmatchup in Matchup_object_list:
                    m_layer = floatmatchup.subset(layer)
                    if m_layer.number() > 0:
                        REF_LAYER_MEAN.append(m_layer.Ref.mean())
                        MODEL_LAYER_MEAN.append(m_layer.Model.mean())
    
                NPOINTS[ivar, iFrame, isub, ilayer] = len(MODEL_LAYER_MEAN)
                if len(MODEL_LAYER_MEAN) > 0:
                    M_LAYER = matchup(np.array(MODEL_LAYER_MEAN), np.array(REF_LAYER_MEAN))
                    BIAS[ivar, iFrame, isub, ilayer] = M_LAYER.bias()
                    RMSE[ivar, iFrame, isub, ilayer] = M_LAYER.RMSE()
                    MOD[ivar, iFrame, isub, ilayer] = M_LAYER.Model.mean()
                    REF[ivar, iFrame, isub, ilayer] = M_LAYER.Ref.mean()



ncOUT = netCDF4.Dataset(args.outfile,'w') #manca l'array times
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
ncvar=ncOUT.createVariable('model', 'f', ('var','time', 'sub','depth'))
ncvar[:] = MOD
ncvar=ncOUT.createVariable('ref', 'f', ('var','time', 'sub','depth'))
ncvar[:] = REF

ncOUT.close()
