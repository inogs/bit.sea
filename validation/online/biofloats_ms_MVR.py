import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates float validation data on a single week, from one single chain run.
    Week is centered on thursday.
    Produces a single file, containing bias, rmse and number of measurements for each subbasin and layer
    for chlorophyll, nitrate and oxygen.
    In this approach we define the measurement as mean on layer of the float profile values.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--date','-d',
                                type = str,
                                required = True,
                                help = 'Date in yyyymmdd format, corresponding to the tuesday = rundate - 7 days')

    parser.add_argument(   '--basedir','-b',
                                type = str,
                                required = True,
                                help = 'dir containing PROFILES/')

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

    parser.add_argument(   '--basedir_clim', '-c',
                                type = str,
                                default = None,
                                required = False,
                                help = "")


    return parser.parse_args()

args = argument()

from commons import timerequestors
from basins import V2 as OGS
from instruments import superfloat as bio_float
from instruments.var_conversions import FLOATVARS
from instruments.matchup_manager import Matchup_Manager
from instruments import check
from commons.Timelist import TimeList
from commons.mask import Mask
from commons.layer import Layer
import numpy as np
from matchup.statistics import matchup
import datetime
import scipy.io.netcdf as NC
from commons.utils import addsep
from basins.region import Rectangle

TheMask  = Mask(args.maskfile)
BASEDIR =  addsep(args.basedir)
if args.basedir_clim is not None:
    BASEDIR_CLIM =  addsep(args.basedir_clim)
outfile  = args.outfile
datestr  = args.date
d = datetime.datetime.strptime(datestr,'%Y%m%d')
R=timerequestors.Weekly_req(d.year,d.month,d.day)

Check_Obj = check.check("",verboselevel=0)

LAYERLIST=[Layer(0,10), Layer(10,30), Layer(30,60), Layer(60,100), Layer(100,150), Layer(150,300)]
VARLIST = ['P_l']
nSub   = len(OGS.MVR.basin_list)
nDepth = len(LAYERLIST)
nVar   = len(VARLIST)

NPOINTS = np.zeros((nVar,nSub,nDepth), np.int32)
NPROFILES = np.zeros((nVar,nSub,nDepth), np.int32)

MODEL_MEAN     = np.zeros((nVar,nSub,nDepth), np.float32)
MODEL_MEAN[:]  = np.nan

REF___MEAN     = np.zeros((nVar,nSub,nDepth), np.float32)
REF___MEAN[:]  = np.nan

MODEL_VARIANCE = np.zeros((nVar,nSub,nDepth), np.float32)
MODEL_VARIANCE[:] = np.nan

REF___VARIANCE = np.zeros((nVar,nSub,nDepth), np.float32)
REF___VARIANCE[:] = np.nan

BIAS           = np.zeros((nVar,nSub,nDepth), np.float32)
RMSE           = np.zeros((nVar,nSub,nDepth), np.float32)
CORR           = np.zeros((nVar,nSub,nDepth), np.float32)
ANOMALY_CORR   = np.zeros((nVar,nSub,nDepth), np.float32)
BIAS[:]        = np.nan
RMSE[:]        = np.nan
CORR[:]        = np.nan
ANOMALY_CORR[:]= np.nan

TI =R.time_interval
TL = TimeList.fromfilenames(TI, BASEDIR + "/PROFILES", "ave*nc")

ALL_PROFILES = bio_float.FloatSelector(None,TI, Rectangle(-6,36,30,46))
M = Matchup_Manager(ALL_PROFILES,TL,BASEDIR)

if args.basedir_clim is not None:
    TLclim = TimeList.fromfilenames(None, BASEDIR_CLIM + "/PROFILES", "ave*nc")
    TLclim.inputFrequency = 'monthly'
    Mclim = Matchup_Manager(ALL_PROFILES,TLclim,BASEDIR_CLIM)


for ivar, var in enumerate(VARLIST):
    print var
    for isub, sub in enumerate(OGS.MVR):
        Profilelist_all = bio_float.FloatSelector(FLOATVARS[var], TI, sub)
        Profilelist = bio_float.remove_bad_sensors(Profilelist_all,FLOATVARS[var])

        Matchup_object_list=[]
        Matchup_object_list_clim=[]
        for p in Profilelist:
            floatmatchup =  M.getMatchups2([p], TheMask.zlevels, var, checkobj=Check_Obj, interpolation_on_Float=True)
            Matchup_object_list.append(floatmatchup)

            if args.basedir_clim is not None:
                floatmatchup_clim =  Mclim.getMatchups2([p], TheMask.zlevels, var,checkobj=Check_Obj, interpolation_on_Float=True)
                Matchup_object_list_clim.append(floatmatchup_clim)
            else:
                Matchup_object_list_clim.append(floatmatchup)

        for ilayer, layer in enumerate(LAYERLIST):
            MODEL_LAYER_MEAN = [] # one value for each suitable profile in (subbasin, layer)
            REF_LAYER_MEAN   = []
            CLIMA_LAYER_MEAN = []
            for floatmatchup,floatmatchup_clim in zip(Matchup_object_list,Matchup_object_list_clim):
                m_layer = floatmatchup.subset(layer)
                if args.basedir_clim is not None:
                    m_layer_clim = floatmatchup_clim.subset(layer)

                if m_layer.number() > 0:
                    REF_LAYER_MEAN.append(m_layer.Ref.mean())
                    MODEL_LAYER_MEAN.append(m_layer.Model.mean())
                    if args.basedir_clim is not None:
                        CLIMA_LAYER_MEAN.append(m_layer_clim.Model.mean())
                    NPOINTS[ivar,isub,ilayer] += len(m_layer.Depth)

            NPROFILES[ivar,isub,ilayer] = len(MODEL_LAYER_MEAN)
            if len(MODEL_LAYER_MEAN) > 0:
                M_LAYER = matchup(np.array(MODEL_LAYER_MEAN), np.array(REF_LAYER_MEAN))
                BIAS[ivar,isub,ilayer] = M_LAYER.bias()
                RMSE[ivar,isub,ilayer] = M_LAYER.RMSE()
                MODEL_MEAN[ivar,isub,ilayer] = np.nanmean(MODEL_LAYER_MEAN)
                REF___MEAN[ivar,isub,ilayer] = np.nanmean(REF_LAYER_MEAN)

            if len(MODEL_LAYER_MEAN) > 1:
                vref, vmodel = M_LAYER.variances()
                MODEL_VARIANCE[ivar,isub,ilayer] = vmodel
                REF___VARIANCE[ivar,isub,ilayer] = vref

            if len(MODEL_LAYER_MEAN) > 2:
                CORR[ivar,isub,ilayer] = M_LAYER.correlation()
                if args.basedir_clim is not None:
                    M_LAYER_ANOMALY = matchup(np.array(MODEL_LAYER_MEAN)-np.array(CLIMA_LAYER_MEAN), \
                                          np.array(REF_LAYER_MEAN)-np.array(CLIMA_LAYER_MEAN))
                    ANOMALY_CORR[ivar,isub,ilayer] = M_LAYER_ANOMALY.correlation()

ncOUT = NC.netcdf_file(outfile,'w')

ncOUT.createDimension('var', nVar)
ncOUT.createDimension('sub', nSub)
ncOUT.createDimension('depth',nDepth)
s=''
for var in VARLIST: s= s+var + ","
setattr(ncOUT, 'varlist',s[:-1])
s=''
for sub in OGS.MVR: s =s+sub.name + ","
setattr(ncOUT,'sublist',s[:-1])
s=''
for layer in LAYERLIST: s =s+layer.string() + ","
setattr(ncOUT,'layerlist',s[:-1])

ncvar=ncOUT.createVariable('bias', 'f', ('var','sub','depth'))
ncvar[:] = BIAS
ncvar=ncOUT.createVariable('rmse', 'f', ('var','sub','depth'))
ncvar[:] = RMSE
ncvar=ncOUT.createVariable('npoints', 'i', ('var','sub','depth'))
ncvar[:] = NPROFILES
ncvar=ncOUT.createVariable('nobs', 'i', ('var','sub','depth'))
ncvar[:] = NPOINTS
ncvar=ncOUT.createVariable('modelmean', 'f', ('var','sub','depth'))
ncvar[:] = MODEL_MEAN
ncvar=ncOUT.createVariable('refmean', 'f', ('var','sub','depth'))
ncvar[:] = REF___MEAN
ncvar=ncOUT.createVariable('modelvar', 'f', ('var','sub','depth'))
ncvar[:] = MODEL_VARIANCE
ncvar=ncOUT.createVariable('refvar', 'f', ('var','sub','depth'))
ncvar[:] = REF___VARIANCE
ncvar=ncOUT.createVariable('corr', 'f', ('var','sub','depth'))
ncvar[:] = CORR
if args.basedir_clim is not None:
    ncvar=ncOUT.createVariable('anomaly_corr', 'f', ('var','sub','depth'))
    ncvar[:] = ANOMALY_CORR

ncOUT.close()
