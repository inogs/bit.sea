import numpy as np
import argparse
import pickle
import matplotlib.pyplot as plt
import scipy.io.netcdf as NC
from datetime import datetime, timedelta
#from Sat import SatManager as Sat
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from commons.dataextractor import DataExtractor
from commons.mask import Mask
from commons.submask import SubMask
from basins import V2
#from postproc import masks
from commons.utils import addsep
#import os


def argument():
    parser = argparse.ArgumentParser(description = '''
    Plot model and sat climatology
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                help = ''' PKL  directory'''
                                )

    parser.add_argument(   '--satclimdir', '-c',
                                type = str,
                                required = True,
                                help = ''' Climatology statistics directory'''
                                )

    parser.add_argument(   '--insat', '-s',
                                type = str,
                                required = True,
                                help = ''' Satellite directory'''
                                )


    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help = ''' File of the model mask
                                           '''
                                )

    parser.add_argument(   '--fileout', '-f',
                                type = str,
                                required = True,
                                help = ''' File of the model mask
                                           '''
                                )

    return parser.parse_args()


args = argument()


CLIMDIR = addsep(args.satclimdir)
INSAT   = addsep(args.insat)
OUTDIR  = addsep(args.outdir)

TheMask = Mask(args.maskfile)
mask200 = TheMask.mask_at_level(200)

print 'Subbasin masks'
submask = {}
subm_open = {}
subnames = ''
for sub in V2.Pred:
    print sub.name
    submask[sub.name] = SubMask(sub,maskobject=TheMask).mask[0,:,:]
    subm_open[sub.name] = submask[sub.name] | (mask200==True)
    subnames += sub.name + ', '

subnames += 'med, '
submask['med'] = np.zeros_like(submask[sub.name])
for sub in V2.Pred:
    submask['med'] = submask['med'] | submask[sub.name]
subm_open['med'] = submask['med'] | (mask200==True)



Timestart = "19000101"
#Timestart = "20160101"
Time__end = "20500101"

TI = TimeInterval(Timestart,Time__end,"%Y%m%d")
TLsat = TimeList.fromfilenames(TI, INSAT, \
            "*.nc", prefix='',dateformat='%Y%m')


Ndates = TLsat.nTimes


# Check on sublist of clima and sat statistics
filestats0 = CLIMDIR + 'statsday001_clim.nc'
NCclim = NC.netcdf_file(filestats0,'r')
# filesat0 = TLsat.filelist[0]
# NCsat = NC.netcdf_file(filestats0,'r')
if (not(NCclim.sublist==subnames)): # or (not(NCsat.sublist==subnames)):
    print '''List of subbasin used for climatology statistics 
             not equal to the actual one'''
    print 'CLIMA list:  ' + NCclim.sublist
    # print 'SAT list:    ' + NCsat.sublist
    print 'ACTUAL list: ' + subnames
    print ' ---  EXIT!  ---'
    NCclim.close()
    import sys
    sys.exit()

NCclim.close()


plotDates = []
satMean = {}
#modelStd = {}
climaMean = {}
climaStd = {}
for sub in V2.P:
    satMean[sub.name] = np.zeros((2,Ndates))
    #modelStd[sub.name] = np.zeros(Ndates)
    climaMean[sub.name] = np.zeros((2,Ndates))
    climaStd[sub.name] = np.zeros((2,Ndates))

for ifile,filesat in enumerate(TLsat.filelist):
    datefile = TLsat.Timelist[ifile]
    plotDates.append(datefile)
    print ' Reading sat file ' + np.str(ifile+1) + ' of ' + np.str(Ndates)
    print '   ... date ' + datetime.strftime(datefile,'%Y-%m-%d')
    De = DataExtractor(TheMask,filename=filesat,varname='lchlm',dimvar=2)
    surfchl = De.filled_values
    
    julianday = int(datefile.strftime('%j'))
    if julianday==366:
        julianday = 365
    fileclima = CLIMDIR + 'statsday' + '%03d' %julianday + '_climAll.nc'
    NCclim = NC.netcdf_file(fileclima,'r')
    climastats = NCclim.variables['SUBstatistics'].data.copy()
    NCclim.close()

    fileclima = CLIMDIR + 'statsday' + '%03d' %julianday + '_climOpen.nc'
    NCclim = NC.netcdf_file(fileclima,'r')
    climast_open = NCclim.variables['SUBstatistics'].data.copy()
    NCclim.close()

    for isub,sub in enumerate(V2.P):
        satMean[sub.name][0,ifile] = np.nanmean(surfchl[submask[sub.name]])
        satMean[sub.name][1,ifile] = np.nanmean(surfchl[subm_open[sub.name]])
        #satStd[sub.name][ifile] = np.nanstd(surfchl[submask[sub.name]])
        climaMean[sub.name][0,ifile] = climastats[isub,0]
        climaMean[sub.name][1,ifile] = climast_open[isub,0]
        climaStd[sub.name][0,ifile] = climastats[isub,9] #average std in the subbasin  
        climaStd[sub.name][1,ifile] = climast_open[isub,9] #average std in the subbasin  


LIST = [i for i in range(4)]

LIST[0] = plotDates
LIST[1] = satMean
LIST[2] = climaMean
LIST[3] = climaStd


filestats = OUTDIR + '/' + args.fileout +'.pkl'


fid = open(filestats,'wb')
pickle.dump(LIST,fid)
fid.close()






