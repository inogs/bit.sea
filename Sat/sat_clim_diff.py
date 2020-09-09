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

    parser.add_argument(   '--climfile', '-c',
                                type = str,
                                required = True,
                                help = ''' Climatology .nc file used to apply check on sat data, e.g. /gss/gss_work/DRES_OGS_BiGe/Observations/CLIMATOLOGY/SAT/CCI02/SatClimatology.nc'''

                                )

    parser.add_argument(   '--starttime', '-s',
                                type = str,
                                required = True,
                                help = ''' Start time ''' 

                                )

    parser.add_argument(   '--endtime', '-e',
                                type = str,
                                required = True,
                                help = ''' End time ''' 

                                )

    parser.add_argument(   '--mesh', '-m',
                                type = str,
                                required = True,
                                choices = ['SatOrigMesh','V4mesh','V1mesh','KD490mesh','SAT1km_mesh', 'Mesh24'],
                                help = ''' Name of the mesh of sat ORIG and used to dump checked data.'''
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
CLIM_FILE = args.climfile

Timestart = args.starttime
Time__end = args.endtime

maskSat = getattr(masks,args.mesh)


TI = TimeInterval(Timestart,Time__end,"%Y%m%d")
TL_orig = TimeList.fromfilenames(TI, ORIGDIR ,"*.nc",prefix='',dateformat='%Y%m%d')

LISTnump = []
print('Read climatology')
MEAN,STD = Sat.readClimatology(CLIM_FILE)

for iTime,filename in enumerate(TL_orig.filelist):
    outfile = OUTDIR + '/DIFF/' + os.path.basename(filename)
    outfilen = OUTDIR + '/DIFF_N/' + os.path.basename(filename)
    iDate = TL_orig.Timelist[iTime] 
    date8 = iDate.strftime('%Y%m%d')
    print date8


    julian = int( iDate.strftime("%j") )
    print(' ... day ' + np.str(julian) + '  of ' + np.str(iDate.year))
    if julian == 366:
        julian = 365


    DAILY_REF_MEAN = MEAN[julian-1,:,:]
    DAILY_REF_STD  =  STD[julian-1,:,:]    
    
    CHL_IN = Sat.readfromfile(filename)

    if args.mesh == 'SatOrigMesh': CHL_IN[581:,164:] = Sat.fillValue # BLACK SEA

    cloudsLandTIME = CHL_IN         == Sat.fillValue
    cloudlandsCLIM = DAILY_REF_MEAN == Sat.fillValue
    LISTnump.append(np.sum(cloudsLandTIME==False))

    CHL_OUT = np.zeros_like(CHL_IN)
    CHL_OUT = CHL_IN-DAILY_REF_MEAN
    CHL_OUT[cloudsLandTIME] = Sat.fillValue
    CHL_OUT[cloudlandsCLIM] = Sat.fillValue
    
    Sat.dumpGenericNativefile(outfile, CHL_OUT, "CHL",mesh=maskSat)

    CHL_OUT_N = np.zeros_like(CHL_IN)
    CHL_OUT_N = (CHL_IN-DAILY_REF_MEAN)/DAILY_REF_STD
    CHL_OUT_N[cloudsLandTIME] = Sat.fillValue
    CHL_OUT_N[cloudlandsCLIM] = Sat.fillValue
    

    Sat.dumpGenericNativefile(outfilen, CHL_OUT_N, "CHL",mesh=maskSat)

filenump = OUTDIR + '/nump' + Timestart + '_' + Time__end + '.txt'

np.savetxt(filenump,LISTnump)
