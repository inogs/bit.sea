import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates validation data on a single day, from one single chain run.
    It works both with weekly and daily satellite data.
    ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument(   '--modelfile','-f',
                                type = str,
                                required = True,
                                help = 'input model file')

    parser.add_argument(   '--satfile','-s',
                                type = str,
                                required = True)
    
    parser.add_argument(   '--climafile','-c',
                                type = str,
                                default = None,
                                required = False)

    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                default = "/marconi/home/usera07ogs/a07ogs00/OPA/V3C/etc/static-data/MED24_125/meshmask.nc",
                                required = False,
                                help = ''' Path of maskfile''')

    parser.add_argument(   '--outfile', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "")

    return parser.parse_args()

args = argument()



import numpy as np
import Sat.SatManager as Sat
import matchup.matchup as matchup
from commons.dataextractor import DataExtractor
from commons.mask import Mask
from commons.submask import SubMask
from basins import V2 as OGS
import glob
import scipy.io.netcdf as NC

def weighted_mean(Conc, Weight):

    Weight_sum      = Weight.sum()
    Mass            = (Conc * Weight).sum()
    Weighted_Mean   = Mass/Weight_sum
    return Weighted_Mean

def weighted_var(Conc, Weight):
    ConcMean = weighted_mean(Conc,Weight)
    Weight_sum      = Weight.sum()
    Mass            = ((Conc-ConcMean)**2 * Weight).sum()
    Weighted_Var    = Mass/Weight_sum
    return Weighted_Var



modfile =args.modelfile
if args.climafile is not None:
    climafile = args.climafile
satfile=args.satfile
outfile = args.outfile

TheMask=Mask(args.maskfile)


nSUB = len(OGS.P.basin_list)

jpk,jpj,jpi =TheMask.shape
mask200_2D = TheMask.mask_at_level(200.0)
dtype = [(sub.name, np.bool) for sub in OGS.P]
SUB = np.zeros((jpj,jpi),dtype=dtype)

for sub in OGS.Pred:
    SUB[sub.name]  = SubMask(sub,maskobject=TheMask).mask_at_level(0)
    if 'atl' in sub.name: continue
    SUB['med'] = SUB['med'] | SUB[sub.name]



COASTNESS_LIST=['coast','open_sea','everywhere']
nCOAST = len(COASTNESS_LIST)
dtype = [(coast, np.bool) for coast in COASTNESS_LIST]
COASTNESS = np.ones((jpj,jpi),dtype=dtype)
COASTNESS['coast']     = ~mask200_2D
COASTNESS['open_sea']  =  mask200_2D


De = DataExtractor(TheMask,filename=modfile,varname='P_l',dimvar=2)
Model = De.values[:,:]

if args.climafile is not None:
    Declim = DataExtractor(TheMask,filename=climafile,varname='P_l',dimvar=2)
    Clima = Declim.values[:,:]
try:
    Sat24 = Sat.readfromfile(satfile,'CHL') # weekly
except:
    Sat24 = Sat.convertinV4format(Sat.readfromfile(satfile, 'CHL'))  # daily

cloudsLand = np.isnan(Sat24)
Sat24[cloudsLand] = -999.0
cloudsLand = Sat24==-999.0
modelLand  = Model > 1.0e+19
nodata     = cloudsLand | modelLand



VALID_POINTS                       = np.zeros((nSUB,nCOAST),np.float32)

MODEL_MEAN                         = np.zeros((nSUB,nCOAST),np.float32)
SAT___MEAN                         = np.zeros((nSUB,nCOAST),np.float32)
MODEL_VARIANCE                     = np.zeros((nSUB,nCOAST),np.float32)
SAT___VARIANCE                     = np.zeros((nSUB,nCOAST),np.float32)
BGC_CLASS4_CHL_BIAS_SURF_BASIN     = np.zeros((nSUB,nCOAST),np.float32)
BGC_CLASS4_CHL_RMS_SURF_BASIN      = np.zeros((nSUB,nCOAST),np.float32)
BGC_CLASS4_CHL_CORR_SURF_BASIN     = np.zeros((nSUB,nCOAST),np.float32)

MODEL_MEAN_LOG                     = np.zeros((nSUB,nCOAST),np.float32)
SAT___MEAN_LOG                     = np.zeros((nSUB,nCOAST),np.float32)
MODEL_VARIANCE_LOG                 = np.zeros((nSUB,nCOAST),np.float32)
SAT___VARIANCE_LOG                 = np.zeros((nSUB,nCOAST),np.float32)
BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG  = np.zeros((nSUB,nCOAST),np.float32)
BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG = np.zeros((nSUB,nCOAST),np.float32)
BGC_CLASS4_CHL_CORR_SURF_BASIN_LOG = np.zeros((nSUB,nCOAST),np.float32)
# varianceCheck = np.zeros((nSUB,nCOAST),np.float32)

ANOMALY_CORR     = np.zeros((nSUB,nCOAST),np.float32)
ANOMALY_CORR_LOG = np.zeros((nSUB,nCOAST),np.float32)

for icoast, coast in enumerate(COASTNESS_LIST):

    for isub, sub in enumerate(OGS.P):
        selection = SUB[sub.name] & (~nodata) & COASTNESS[coast]
        M = matchup.matchup(Model[selection], Sat24[selection])
        VALID_POINTS[isub,icoast] = M.number()
        if M.number() > 0 :
            BGC_CLASS4_CHL_RMS_SURF_BASIN[isub,icoast]  = M.RMSE()
            BGC_CLASS4_CHL_BIAS_SURF_BASIN[isub,icoast] = M.bias()
            BGC_CLASS4_CHL_CORR_SURF_BASIN[isub,icoast] = M.correlation()

            weight = TheMask.area[selection]
            MODEL_MEAN[isub,icoast] = weighted_mean( M.Model,weight)
            SAT___MEAN[isub,icoast] = weighted_mean( M.Ref,  weight)
            MODEL_VARIANCE[isub,icoast] = weighted_var( M.Model,weight)
            SAT___VARIANCE[isub,icoast] = weighted_var( M.Ref,  weight)

            Mlog = matchup.matchup(np.log10(Model[selection]), np.log10(Sat24[selection])) #add matchup based on logarithm
            BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG[isub,icoast]  = Mlog.RMSE()
            BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG[isub,icoast] = Mlog.bias()
            BGC_CLASS4_CHL_CORR_SURF_BASIN_LOG[isub,icoast] = Mlog.correlation()

            MODEL_MEAN_LOG[isub,icoast] = weighted_mean( Mlog.Model,weight)
            SAT___MEAN_LOG[isub,icoast] = weighted_mean( Mlog.Ref,  weight)
            MODEL_VARIANCE_LOG[isub,icoast] = weighted_var( Mlog.Model,weight)
            SAT___VARIANCE_LOG[isub,icoast] = weighted_var( Mlog.Ref,  weight)

            if args.climafile is not None:
                Mclima = matchup.matchup(Model[selection]-Clima[selection],  \
                                        Sat24[selection]-Clima[selection])
                ANOMALY_CORR[isub,icoast] = Mclima.correlation()

                Mclima_log = matchup.matchup(np.log10(Model[selection])-np.log10(Clima[selection]), \
                                            np.log10(Sat24[selection])-np.log10(Clima[selection]))
                ANOMALY_CORR_LOG[isub,icoast] = Mclima_log.correlation()
        else:
            BGC_CLASS4_CHL_RMS_SURF_BASIN[isub,icoast] = np.nan
            BGC_CLASS4_CHL_BIAS_SURF_BASIN[isub,icoast]= np.nan
            BGC_CLASS4_CHL_CORR_SURF_BASIN[isub,icoast]= np.nan
            MODEL_MEAN[isub,icoast] = np.nan
            SAT___MEAN[isub,icoast] = np.nan
            MODEL_VARIANCE[isub,icoast] = np.nan
            SAT___VARIANCE[isub,icoast] = np.nan
            
            BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG[isub,icoast] = np.nan
            BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG[isub,icoast]= np.nan
            BGC_CLASS4_CHL_CORR_SURF_BASIN_LOG[isub,icoast]= np.nan
            MODEL_MEAN_LOG[isub,icoast] = np.nan
            SAT___MEAN_LOG[isub,icoast] = np.nan
            MODEL_VARIANCE_LOG[isub,icoast] = np.nan
            SAT___VARIANCE_LOG[isub,icoast] = np.nan

            ANOMALY_CORR[isub,icoast]     = np.nan
            ANOMALY_CORR_LOG[isub,icoast] = np.nan


ncOUT = NC.netcdf_file(outfile,'w')
ncOUT.createDimension('nsub', nSUB)
ncOUT.createDimension('ncoast', nCOAST)
s=''
for sub in OGS.P: s =s+sub.name + ","
setattr(ncOUT,'sublist',s)
s=''
for coast in COASTNESS_LIST: s =s+coast + ","
setattr(ncOUT,'coastlist',s)


ncvar= ncOUT.createVariable('number', 'i', ('nsub','ncoast'))
ncvar[:]= VALID_POINTS

ncvar = ncOUT.createVariable('BGC_CLASS4_CHL_RMS_SURF_BASIN', 'f', ('nsub','ncoast'))
ncvar[:] = BGC_CLASS4_CHL_RMS_SURF_BASIN

ncvar = ncOUT.createVariable('BGC_CLASS4_CHL_BIAS_SURF_BASIN', 'f', ('nsub','ncoast'))
ncvar[:] = BGC_CLASS4_CHL_BIAS_SURF_BASIN

ncvar = ncOUT.createVariable('BGC_CLASS4_CHL_CORR_SURF_BASIN', 'f', ('nsub','ncoast'))
ncvar[:] = BGC_CLASS4_CHL_CORR_SURF_BASIN

ncvar = ncOUT.createVariable('MODEL_MEAN', 'f', ('nsub','ncoast'))
ncvar[:] = MODEL_MEAN

ncvar=ncOUT.createVariable('SAT___MEAN', 'f', ('nsub','ncoast'))
ncvar[:] = SAT___MEAN

ncvar = ncOUT.createVariable('MODEL_VARIANCE', 'f', ('nsub','ncoast'))
ncvar[:] = MODEL_VARIANCE

ncvar=ncOUT.createVariable('SAT___VARIANCE', 'f', ('nsub','ncoast'))
ncvar[:] = SAT___VARIANCE

ncvar = ncOUT.createVariable('BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG', 'f', ('nsub','ncoast'))
ncvar[:] = BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG

ncvar= ncOUT.createVariable('BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG', 'f', ('nsub','ncoast'))
ncvar[:] =BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG

ncvar= ncOUT.createVariable('BGC_CLASS4_CHL_CORR_SURF_BASIN_LOG', 'f', ('nsub','ncoast'))
ncvar[:] =BGC_CLASS4_CHL_CORR_SURF_BASIN_LOG

ncvar = ncOUT.createVariable('MODEL_MEAN_LOG', 'f', ('nsub','ncoast'))
ncvar[:] = MODEL_MEAN_LOG

ncvar=ncOUT.createVariable('SAT___MEAN_LOG', 'f', ('nsub','ncoast'))
ncvar[:] = SAT___MEAN_LOG

ncvar = ncOUT.createVariable('MODEL_VARIANCE_LOG', 'f', ('nsub','ncoast'))
ncvar[:] = MODEL_VARIANCE_LOG

ncvar=ncOUT.createVariable('SAT___VARIANCE_LOG', 'f', ('nsub','ncoast'))
ncvar[:] = SAT___VARIANCE_LOG

if args.climafile is not None:
    ncvar=ncOUT.createVariable('ANOMALY_CORRELATION', 'f', ('nsub','ncoast'))
    ncvar[:] = ANOMALY_CORR

    ncvar=ncOUT.createVariable('ANOMALY_CORRELATION_LOG', 'f', ('nsub','ncoast'))
    ncvar[:] = ANOMALY_CORR_LOG

ncOUT.close()
    

    
    
