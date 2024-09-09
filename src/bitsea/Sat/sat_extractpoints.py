import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Calculates differences wrt clima 
    also normalized eith clima std
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--origdir', '-i',
                                type = str,
                                required = True,
                                help = ''' ORIG sat directory, e.g. /gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MODIS/DAILY/ORIG/'''
                                )

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                help = ''' outdir with diff files ''' 
                                )

#    parser.add_argument(   '--starttime', '-s',
#                                type = str,
#                                required = True,
#                                help = ''' Start time ''' 
#                                )

#    parser.add_argument(   '--endtime', '-e',
#                                type = str,
#                                required = True,
#                                help = ''' End time ''' 
#                                )

    parser.add_argument(   '--mesh', '-m',
                                type = str,
                                required = True,
                                choices = ['SatOrigMesh','V4mesh','V1mesh','KD490mesh','SAT1km_mesh', 'Mesh24'],
                                help = ''' Name of the mesh of sat ORIG and used to dump checked data.'''
                                )

    parser.add_argument(   '--lon', '-n',
                                type = float,
                                required = True,
                                nargs='+',
                                help = ''' Longitudes, same number of lats than lons''' 
                                )

    parser.add_argument(   '--lat', '-t',
                                type = float,
                                required = True,
                                nargs='+',
                                help = ''' Latitudes, same number of lats than lons''' 

                                )
    return parser.parse_args()


args = argument()
import numpy as np
from commons.time_interval import TimeInterval
from commons.Timelist import TimeList
from commons.utils import addsep
from postproc import masks
import os

import SatManager as Sat


ORIGDIR   = addsep(args.origdir)
OUTDIR  = addsep(args.outdir)

#Timestart = args.starttime
#Time__end = args.endtime
Timestart = '19500101'
Time__end = '20500101'

maskSat = getattr(masks,args.mesh)
lonLIST = args.lon
latLIST = args.lat

ilonLIST = []
ilatLIST = []
labelLIST = []
for ii,lonp in enumerate(lonLIST):
    ilonp = np.argmin(abs(lonp-maskSat.lon))
    ilonLIST.append(ilonp)

    latp = latLIST[ii]
    ilatp = np.argmin(abs(latp-maskSat.lat))
    ilatLIST.append(ilatp)

    label = str(lonp) + '_' + str(latp)
    labelLIST.append(label)

TI = TimeInterval(Timestart,Time__end,"%Y%m%d")
TL_orig = TimeList.fromfilenames(TI, ORIGDIR ,"*.nc",prefix='',dateformat='%Y%m%d')
#TL_orig = TimeList.fromfilenames(None, ORIGDIR ,"*.nc",prefix='',dateformat='%Y%m%d')

LISTnump = []

DICTchl = {}
for label in labelLIST:
    DICTchl[label] = []

dates = []

for iTime,filename in enumerate(TL_orig.filelist):
    iDate = TL_orig.Timelist[iTime] 
    date8 = iDate.strftime('%Y%m%d')
    print date8
    dates.append(date8)
    
    if args.mesh=='V4mesh':
        CHL_IN = Sat.readfromfile(filename,var='lchlm')
    else:
        CHL_IN = Sat.readfromfile(filename)
    for iip,ilonp in enumerate(ilonLIST):
        ilatp = ilatLIST[iip]
        chlp = CHL_IN[ilatp,ilonp]
        if chlp>0:
            DICTchl[labelLIST[iip]].append(chlp)
        else:
            DICTchl[labelLIST[iip]].append(np.nan)
        


for label in labelLIST:
    filep = OUTDIR + '/satp' + label  + '.txt'
    np.savetxt(filep,DICTchl[label])


filed = OUTDIR + '/dates.txt'
np.savetxt(filed,dates,fmt='%s')

