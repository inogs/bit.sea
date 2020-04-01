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



modfile =args.modelfile
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
    SUB['med'] = SUB['med'] | SUB[sub.name]



COASTNESS_LIST=['coast','open_sea','everywhere']
nCOAST = len(COASTNESS_LIST)
dtype = [(coast, np.bool) for coast in COASTNESS_LIST]
COASTNESS = np.ones((jpj,jpi),dtype=dtype)
COASTNESS['coast']     = ~mask200_2D;
COASTNESS['open_sea']  =  mask200_2D;
#COASTNESS['everywhere'] = True;

De = DataExtractor(TheMask,filename=modfile,varname='P_l')
Model = De.values[0,:,:]
try:
#    Sat24 = Sat.readfromfile(satfile,'lchlm') # weekly
    Sat24 = Sat.readfromfile(satfile,'CHL') # weekly
except:
    Sat24 = Sat.convertinV4format(Sat.readfromfile(satfile, 'CHL'))  # daily

cloudsLand = np.isnan(Sat24)
Sat24[cloudsLand] = -999.0
cloudsLand = Sat24==-999.0; 
modelLand  = Model > 1.0e+19
nodata     = cloudsLand | modelLand



BGC_CLASS4_CHL_RMS_SURF_BASIN      = np.zeros((nSUB,nCOAST),np.float32)
BGC_CLASS4_CHL_BIAS_SURF_BASIN     = np.zeros((nSUB,nCOAST),np.float32)
VALID_POINTS                       = np.zeros((nSUB,nCOAST),np.float32)
MODEL_MEAN                         = np.zeros((nSUB,nCOAST),np.float32)
SAT___MEAN                         = np.zeros((nSUB,nCOAST),np.float32)
BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG  = np.zeros((nSUB,nCOAST),np.float32)
BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG = np.zeros((nSUB,nCOAST),np.float32)

for icoast, coast in enumerate(COASTNESS_LIST):

    for isub, sub in enumerate(OGS.P):
        selection = SUB[sub.name] & (~nodata) & COASTNESS[coast]
        M = matchup.matchup(Model[selection], Sat24[selection])
        VALID_POINTS[isub,icoast] = M.number()
        if M.number() > 0 :
            BGC_CLASS4_CHL_RMS_SURF_BASIN[isub,icoast]  = M.RMSE()
            BGC_CLASS4_CHL_BIAS_SURF_BASIN[isub,icoast] = M.bias()
                 
            weight = TheMask.area[selection]
            MODEL_MEAN[isub,icoast] = weighted_mean( M.Model,weight)
            SAT___MEAN[isub,icoast] = weighted_mean( M.Ref,  weight)
        
            Mlog = matchup.matchup(np.log10(Model[selection]), np.log10(Sat24[selection])) #add matchup based on logarithm
            BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG[isub,icoast]  = Mlog.RMSE()
            BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG[isub,icoast] = Mlog.bias()
        else:
            BGC_CLASS4_CHL_RMS_SURF_BASIN[isub,icoast] = np.nan
            BGC_CLASS4_CHL_BIAS_SURF_BASIN[isub,icoast]= np.nan
            MODEL_MEAN[isub,icoast] = np.nan
            SAT___MEAN[isub,icoast] = np.nan
            BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG[isub,icoast] = np.nan
            BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG[isub,icoast]= np.nan



ncOUT = NC.netcdf_file(outfile,'w')
ncOUT.createDimension('nsub', nSUB)
ncOUT.createDimension('ncoast', nCOAST)
s='';
for sub in OGS.P: s =s+sub.name + ","
setattr(ncOUT,'sublist',s)
s='';
for coast in COASTNESS_LIST: s =s+coast + ","
setattr(ncOUT,'coastlist',s)


ncvar= ncOUT.createVariable('number', 'i', ('nsub','ncoast'))
ncvar[:]= VALID_POINTS

ncvar = ncOUT.createVariable('BGC_CLASS4_CHL_RMS_SURF_BASIN', 'f', ('nsub','ncoast'))
ncvar[:] = BGC_CLASS4_CHL_RMS_SURF_BASIN

ncvar = ncOUT.createVariable('BGC_CLASS4_CHL_BIAS_SURF_BASIN', 'f', ('nsub','ncoast'))
ncvar[:] = BGC_CLASS4_CHL_BIAS_SURF_BASIN

ncvar = ncOUT.createVariable('MODEL_MEAN', 'f', ('nsub','ncoast'))
ncvar[:] = MODEL_MEAN

ncvar=ncOUT.createVariable('SAT___MEAN', 'f', ('nsub','ncoast'))
ncvar[:] = SAT___MEAN

ncvar = ncOUT.createVariable('BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG', 'f', ('nsub','ncoast'))
ncvar[:] = BGC_CLASS4_CHL_RMS_SURF_BASIN_LOG

ncvar= ncOUT.createVariable('BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG', 'f', ('nsub','ncoast'))
ncvar[:] =BGC_CLASS4_CHL_BIAS_SURF_BASIN_LOG

ncOUT.close()
    

    
    
