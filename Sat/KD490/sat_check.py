import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Apply check based on climatology to sat ORIG files for Kd490
    Produces CHECKED files for each date at satellite resolution.
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--origdir', '-i',
                                type = str,
                                required = True,
                                help = ''' ORIG sat directory'''
                                )

    parser.add_argument(   '--checkdir', '-o',
                                type = str,
                                required = True,
                                help = ''' Base for CHECKED sat directory'''
                                )

    parser.add_argument(   '--climfile', '-c',
                                type = str,
                                required = True,
                                help = ''' Climatology .nc file used to apply check on sat data'''
                                )

    parser.add_argument(   '--mesh', '-m',
                                type = str,
                                required = True,
                                choices = ['SatOrigMesh','V4mesh','V1mesh','KD490mesh','SAT1km_mesh', 'Mesh24'],
                                help = ''' Name of the mesh of sat ORIG and used to dump checked data.'''
                                )


    return parser.parse_args()


args = argument()



from commons.Timelist import TimeList
from commons.utils import addsep
from commons.time_interval import TimeInterval
from postproc import masks
import numpy as np
import os

from Sat import SatManager as Sat

#ORIGDIR  ="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/KD490/DAILY/ORIG/"
ORIGDIR = addsep(args.origdir)
#CHECKDIR ="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/KD490/DAILY/CHECKED/"
CHECKDIR = addsep(args.checkdir)
#CLIM_FILE="/gss/gss_work/DRES_OGS_BiGe/Observations/CLIMATOLOGY/SAT/KD490/Climatology_KD490.nc"
CLIM_FILE = args.climfile
maskSat = getattr(masks,args.mesh)

reset = False

Timestart="19990101"
Time__end="20500101"
TI = TimeInterval(Timestart,Time__end,"%Y%m%d")
TL_orig = TimeList.fromfilenames(TI, ORIGDIR ,"*.nc",prefix='',dateformat='%Y%m%d')


somecheck = False
for filename in TL_orig.filelist:
    outfile = CHECKDIR + os.path.basename(filename)
    if not os.path.exists(outfile) :
        somecheck = True
        break

if somecheck:
    MEAN,STD = Sat.readClimatology(CLIM_FILE)
else:
    print "All checks done"

for iTime, filename in enumerate(TL_orig.filelist):
    outfile = CHECKDIR + os.path.basename(filename)
    exit_condition = os.path.exists(outfile) and (not reset) 
    if exit_condition: 
        continue
    julian = int( TL_orig.Timelist[iTime].strftime("%j") )
    if julian==366: julian=365 
    DAILY_REF_MEAN = MEAN[julian-1,:,:]
    DAILY_REF_STD  =  STD[julian-1,:,:]    

    CHL_IN = Sat.readfromfile(filename,'KD490')
    #CHL_IN[581:,164:] = Sat.fillValue # BLACK SEA
    cloudsLandTIME = CHL_IN         == Sat.fillValue
    cloudlandsCLIM = DAILY_REF_MEAN == Sat.fillValue
    
    CHL_OUT = CHL_IN.copy()
    CHL_OUT[cloudsLandTIME] = Sat.fillValue
    CHL_OUT[cloudlandsCLIM] = Sat.fillValue
    counter_refNAN = (~cloudsLandTIME & cloudlandsCLIM).sum(axis=None)
    
    
    outOfRange = np.abs(CHL_IN - DAILY_REF_MEAN) > DAILY_REF_STD *2.0
    outOfRange[cloudsLandTIME | cloudlandsCLIM ] = False
    
    counter_elim = outOfRange.sum(axis = None)
    CHL_OUT[outOfRange] = Sat.fillValue 
    
    print filename
    print 'Rejection:  after check', counter_elim, ' values'
    print 'rejected for NAN in Climatology', counter_refNAN, ' values'
    Sat.dumpGenericNativefile(outfile, CHL_OUT, "KD490",mesh=maskSat)


    
     
    





