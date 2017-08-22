import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Needs a profiler.py, already executed.

    Produces a single file, containing bias, rmse and number of measurements for each subbasin and layer
    for chlorophyll, nitrate and oxygen.
    In this approach we define the measurement as mean on layer of the float profile values.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                default = "/pico/home/usera07ogs/a07ogs00/OPA/V2C/etc/static-data/MED1672_cut/MASK/meshmask.nc",
                                required = False,
                                help = ''' Path of maskfile''')

    parser.add_argument(   '--outfile', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "")

    return parser.parse_args()

args = argument()

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


TheMask  = Mask(args.maskfile)

NO_N3n = ["6903197","6901767","6901773","6901771"]
NO_O2o = ["6901766","6901510"]

LAYERLIST=[Layer(0,10), Layer(10,30), Layer(30,60), Layer(60,100), Layer(100,150), Layer(150,300), Layer(300,600), Layer(600,1000)]
VARLIST = ['P_l','N3n','O2o']
read_adjusted = [True,True,False]
nSub   = len(OGS.NRT3.basin_list)
nDepth = len(LAYERLIST)
nVar   = len(VARLIST)


nFrames = TL.nTimes
BIAS    = np.zeros((nVar,nFrames,nSub,nDepth), np.float32)
RMSE    = np.zeros((nVar,nFrames,nSub,nDepth), np.float32)
NPOINTS = np.zeros((nVar,nFrames, nSub,nDepth), np.int32)

N3n_PROFILES = []
O2o_PROFILES = []

for ip, p in enumerate(ALL_PROFILES):
    if ( p.name() not in NO_N3n): 
	N3n_PROFILES.append(p)
    if ( p.name() not in NO_O2o):
	O2o_PROFILES.append(p)

M_Pl = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)
M_N3n = Matchup_Manager(N3n_PROFILES,TL,BASEDIR)
M_O2o = Matchup_Manager(O2o_PROFILES,TL,BASEDIR)


for iFrame, req in enumerate(TL.getOwnList()):
    if req.time_interval.start_time < TL.timeinterval.start_time : req.time_interval.start_time = TL.timeinterval.start_time
    if req.time_interval.end_time   > TL.timeinterval.end_time   : req.time_interval.end_time   = TL.timeinterval.end_time
    print req.time_interval
    for ivar, var in enumerate(VARLIST):
        print var
	if (var == "P_l"): M = M_Pl
	if (var == "N3n"): M = M_N3n
	if (var == "O2o"): M = M_O2o
        for isub, sub in enumerate(OGS.NRT3):
            Profilelist_ALL = bio_float.FloatSelector(LOVFLOATVARS[var], req.time_interval, sub)
	    Profilelist = []
	    for p in Profilelist_ALL:
		if (  var == "P_l"): Profilelist = Profilelist_ALL
	        if (( var == "N3n" ) and (p.name() not in NO_N3n)): Profilelist.append(p)
                if (( var == "O2o" ) and (p.name() not in NO_O2o)): Profilelist.append(p)
 
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




ncOUT = NC.netcdf_file(args.outfile,'w') #manca l'array times
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
