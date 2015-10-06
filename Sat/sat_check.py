from postproc.Timelist import *
import numpy as np

import SatManager as Sat

ORIGDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MODIS/DAILY/ORIG/"
CHECKDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MODIS/DAILY/CHECKED/"
CLIM_FILE="/pico/home/usera07ogs/a07ogs00/OPA/V4/etc/static-data/DA/CHECKSAT/SatClimatology.nc"

reset = False

Timestart="19500101"
Time__end="20500101"
IonamesFile = '../postproc/IOnames_sat.xml'
TL_orig = TimeList(Timestart,Time__end, ORIGDIR ,"*.nc",IonamesFile)
TLCheck = TimeList(Timestart,Time__end, CHECKDIR,"*.nc",IonamesFile)

ORIG_NAMES=[os.path.basename(i) for i in TL_orig.filelist]
CHECKNAMES=[os.path.basename(i) for i in TLCheck.filelist]

if reset : CHECKNAMES = []


toCheckList=[]
for iTimeorig, filename in enumerate(ORIG_NAMES):
    if filename not in CHECKNAMES:
        toCheckList.append((filename,iTimeorig))

        
if len(toCheckList)>0:
    MEAN,STD = Sat.readClimatology(CLIM_FILE)
else:
    print "All checks done"

for filename, iTimeorig in toCheckList:
    julian = int( TL_orig.Timelist[iTimeorig].strftime("%j") )
    DAILY_REF_MEAN = MEAN[julian-1,:,:]
    DAILY_REF_STD  =  STD[julian-1,:,:]    
    
    CHL_IN = Sat.readfromfile(TL_orig.filelist[iTimeorig])
    CHL_IN[581:,164:] = Sat.fillValue # BLACK SEA
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
    Sat.dumpfile(CHECKDIR + filename, CHL_OUT)


    
     
    





