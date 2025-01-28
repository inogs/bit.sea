import argparse
from bitsea.utilities.argparse_types import existing_dir_path, existing_file_path
from bitsea.utilities.argparse_types import date_from_str, generic_path
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates float validation data on a single week, from one single chain run.
    Week is centered on thursday.
    Produces a single file, containing bias, rmse and number of measurements for each subbasin and layer
    for chlorophyll, nitrate and oxygen.
    In this approach we define the measurement as mean on layer of the float profile values.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--date','-d',
                                type = date_from_str,
                                required = True,
                                help = 'Date in yyyymmdd format, corresponding to the tuesday = rundate - 7 days')

    parser.add_argument(   '--basedir','-b',
                                type = existing_dir_path,
                                required = True,
                                help = """ PROFILATORE dir, already generated""")

    parser.add_argument(   '--maskfile', '-m',
                                type = existing_file_path,
                                required = True)

    parser.add_argument(   '--outfile', '-o',
                                type = generic_path,
                                required = True,
                                help = "")

    return parser.parse_args()

args = argument()

from bitsea.commons import timerequestors
from bitsea.basins import V2 as OGS
from bitsea.instruments import superfloat as bio_float
from bitsea.instruments.var_conversions import FLOATVARS
from bitsea.instruments.matchup_manager import Matchup_Manager
from bitsea.commons.Timelist import TimeList
from bitsea.commons.mask import Mask
from bitsea.commons.layer import Layer
import numpy as np
from bitsea.matchup.statistics import matchup
import datetime
import netCDF4 as NC

from bitsea.basins.region import Rectangle

TheMask  = Mask(args.maskfile)
BASEDIR =  args.basedir
d  = args.date
R=timerequestors.Weekly_req(d.year,d.month,d.day)

LAYERLIST=[Layer(0,10), Layer(10,30), Layer(30,60), Layer(60,100), Layer(100,150), Layer(150,300), Layer(300,600), Layer(600,1000)]
VARLIST = ['P_l','N3n','O2o']
nSub   = len(OGS.NRT3.basin_list)
nDepth = len(LAYERLIST)
nVar   = len(VARLIST)

BIAS    = np.zeros((nVar,nSub,nDepth), np.float32)
RMSE    = np.zeros((nVar,nSub,nDepth), np.float32)
NPOINTS = np.zeros((nVar,nSub,nDepth), np.int32)

Init_time=datetime.datetime(2021,1,1)
TI =R.time_interval
if TI.start_time < Init_time : TI.start_time = Init_time
TL = TimeList.fromfilenames(TI, BASEDIR / "PROFILES", "ave*nc")

ALL_PROFILES = bio_float.FloatSelector(None,TI, Rectangle(-6,36,30,46))
M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)


for ivar, var in enumerate(VARLIST):
    print (var)
    for isub, sub in enumerate(OGS.NRT3):
        Profilelist_all = bio_float.FloatSelector(FLOATVARS[var], TI, sub)
        Profilelist = bio_float.remove_bad_sensors(Profilelist_all,FLOATVARS[var])
        nProfiles = len(Profilelist)

        Matchup_object_list=[]
        for ip in range(nProfiles):
            floatmatchup =  M.getMatchups2([Profilelist[ip]], TheMask.zlevels, var)
            Matchup_object_list.append(floatmatchup)

        for ilayer, layer in enumerate(LAYERLIST):
            MODEL_LAYER_MEAN = [] # one value for each suitable profile in (subbasin, layer)
            REF_LAYER_MEAN   = []
            for floatmatchup in Matchup_object_list:
                m_layer = floatmatchup.subset(layer)

                if m_layer.number() > 0:
                    REF_LAYER_MEAN.append(m_layer.Ref.mean())
                    MODEL_LAYER_MEAN.append(m_layer.Model.mean())

            NPOINTS[ivar,isub,ilayer] = len(MODEL_LAYER_MEAN)
            if len(MODEL_LAYER_MEAN) > 0:
                M_LAYER = matchup(np.array(MODEL_LAYER_MEAN), np.array(REF_LAYER_MEAN))
                BIAS[ivar,isub,ilayer] = M_LAYER.bias()
                RMSE[ivar,isub,ilayer] = M_LAYER.RMSE()



ncOUT = NC.Dataset(args.outfile,'w')

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

ncvar=ncOUT.createVariable('bias', 'f', ('var','sub','depth'))
ncvar[:] = BIAS
ncvar=ncOUT.createVariable('rmse', 'f', ('var','sub','depth'))
ncvar[:] = RMSE
ncvar=ncOUT.createVariable('npoints', 'i', ('var','sub','depth'))
ncvar[:] = NPOINTS

ncOUT.close()
