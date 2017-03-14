from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
import numpy as np
import os

from Sat import SatManager as Sat

ORIGDIR  ="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/KD490/DAILY/ORIG/"
CHECKDIR ="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/KD490/DAILY/CHECKED/"
CLIM_FILE="/gss/gss_work/DRES_OGS_BiGe/Observations/CLIMATOLOGY/SAT/KD490/Climatology_KD490.nc"

reset = False

Timestart="19990101"
Time__end="20160101"
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
    Sat.dump_SAT1km_nativefile(outfile, CHL_OUT)


    
     
    





