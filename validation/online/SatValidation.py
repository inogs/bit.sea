import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates validation data on a single day, from one single chain run.
    It works both with weekly and daily satellite data.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--date','-d',
                                type = str,
                                required = True,
                                help = 'first day of forecast')

    parser.add_argument(   '--forecastdir','-f',
                                type = str,
                                required = True,
                                help = 'input dir validation tmp')
    parser.add_argument(   '--satdir','-s',
                                type = str,
                                required = True,
                                help = 'Satellite directory, usually weekly or daily')

    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                default = "/pico/home/usera07ogs/a07ogs00/OPA/V4/etc/static-data/MED1672_cut/MASK/meshmask.nc",
                                required = False,
                                help = ''' Path of maskfile''')

    parser.add_argument(   '--outfile', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "")

    return parser.parse_args()

args = argument()



from commons.utils import addsep
import numpy as np
import Sat.SatManager as Sat
import matchup.matchup as matchup
from commons import netcdf3
from commons.mask import Mask
from commons.submask import SubMask
from basins import OGS
import glob
import scipy.io.netcdf as NC

def weighted_mean(Conc, Weight):

    Weight_sum      = Weight.sum()
    Mass            = (Conc * Weight).sum()
    Weighted_Mean   = Mass/Weight_sum
    return Weighted_Mean



date=args.date
MODELDIR=addsep(args.forecastdir)
REF_DIR=addsep(args.satdir)
outfile = args.outfile

TheMask=Mask(args.maskfile)


nSUB = len(OGS.P.basin_list)

jpk,jpj,jpi =TheMask.shape
dtype = [(sub.name, np.bool) for sub in OGS.P]
SUB = np.zeros((jpj,jpi),dtype=dtype)
for sub in OGS.P:
    print sub
    sbmask         = SubMask(sub,maskobject=TheMask).mask
    SUB[sub.name]  = sbmask[0,:,:].copy()

mask200_2D = TheMask.mask_at_level(200.0)


satfile = glob.glob(REF_DIR+date+"*")[0]
modfile = MODELDIR + "ave." + date + "-12:00:00.nc" 


Model = netcdf3.read_3d_file(modfile, 'P_l')[0,:,:]
try:
    Sat16 = Sat.readfromfile(satfile,'lchlm') # weekly
except:
    Sat16 = Sat.convertinV4format(Sat.readfromfile(satfile, 'CHL'))  # daily

cloudsLand = np.isnan(Sat16)
Sat16[cloudsLand] = -999.0
cloudsLand = Sat16==-999.0; 
modelLand  = Model > 1.0e+19
nodata     = cloudsLand | modelLand




BGC_CLASS4_CHL_RMS_SURF_BASIN      = np.zeros((nSUB,),np.float32)
BGC_CLASS4_CHL_BIAS_SURF_BASIN     = np.zeros((nSUB,),np.float32)
VALID_POINTS                       = np.zeros((nSUB,),np.float32)
MODEL_MEAN                         = np.zeros((nSUB,),np.float32)
SAT___MEAN                         = np.zeros((nSUB,),np.float32)
BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG  = np.zeros((nSUB,),np.float32)
BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG = np.zeros((nSUB,),np.float32)



for isub, sub in enumerate(OGS.P):
    print sub
    selection = SUB[sub.name] & (~nodata) & mask200_2D
    M = matchup.matchup(Model[selection], Sat16[selection])
    VALID_POINTS[isub] = M.number()
    if M.number() > 0 :
        BGC_CLASS4_CHL_RMS_SURF_BASIN[isub]  = M.RMSE()
        BGC_CLASS4_CHL_BIAS_SURF_BASIN[isub] = M.bias()
             
        weight = TheMask.area[selection]
        MODEL_MEAN[isub] = weighted_mean( M.Model,weight)
        SAT___MEAN[isub] = weighted_mean( M.Ref,  weight)
    
        Mlog = matchup.matchup(np.log10(Model[selection]), np.log10(Sat16[selection])) #add matchup based on logarithm
        BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG[isub]  = Mlog.RMSE()
        BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG[isub] = Mlog.bias()
    else:
        BGC_CLASS4_CHL_RMS_SURF_BASIN = np.nan
        BGC_CLASS4_CHL_BIAS_SURF_BASIN= np.nan
        MODEL_MEAN = np.nan
        SAT___MEAN = np.nan
        BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG = np.nan
        BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG= np.nan



ncOUT = NC.netcdf_file(outfile,'w')
ncOUT.createDimension('nsub', nSUB)
s='';
for sub in OGS.P: s =s+sub.name + ","
setattr(ncOUT,'sublist',s)

ncvar= ncOUT.createVariable('number', 'i', ('nsub',))
ncvar[:]= VALID_POINTS

ncvar = ncOUT.createVariable('BGC_CLASS4_CHL_RMS_SURF_BASIN', 'f', ('nsub',))
ncvar[:] = BGC_CLASS4_CHL_RMS_SURF_BASIN

ncvar = ncOUT.createVariable('BGC_CLASS4_CHL_BIAS_SURF_BASIN', 'f', ('nsub',))
ncvar[:] = BGC_CLASS4_CHL_BIAS_SURF_BASIN

ncvar = ncOUT.createVariable('MODEL_MEAN', 'f', ('nsub',))
ncvar[:] = MODEL_MEAN

ncvar=ncOUT.createVariable('SAT___MEAN', 'f', ('nsub',))
ncvar[:] = SAT___MEAN

ncvar = ncOUT.createVariable('BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG', 'f', ('nsub',))
ncvar[:] = BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG

ncvar= ncOUT.createVariable('BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG', 'f', ('nsub',))
ncvar[:] =BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG

ncOUT.close()
    

    
    