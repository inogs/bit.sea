import postproc
from postproc import Timelist
import glob,os
import numpy as np
import scipy.io.netcdf as NC


ORIGDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MODIS/DAILY/ORIG/"
CHECKDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MODIS/DAILY/CHECKED/"
CLIM_FILE="/pico/home/usera07ogs/a07ogs00/OPA/V4/etc/static-data/DA/CHECKSAT/SatClimatology.nc"


TL_orig = Timelist("19500101","20500101", ORIGDIR ,"*.nc",'postproc/IOnames_sat.xml')
TLCheck = Timelist("19500101","20500101", CHECKDIR,"*.nc",'postproc/IOnames_sat.xml')

#os.chdir(ORIGDIR); ORIG__LIST=glob.glob("*nc")os.chdir(CHECKDIR); CHECKLIST=glob.glob("*nc")

toCheckList=[]

ORIG_NAMES=[os.path.basename(i) for i in TL_orig.filelist]
CHECKNAMES=[os.path.basename(i) for i in TLCheck.filelist]

for it, filename in enumerate(ORIG_NAMES):
    if filename in CHECKNAMES:
        pass
    else:
        toCheckList.append((filename,it))

fillValue=-999.0
        
if len(toCheckList)>0:
    #load Satclimatology
    DAILY_REF_MEAN = 0
    DAILY_REF_STD  = 0. 
    jpi = 300
    jpj = 600
    
    for f in toCheckList:
        filename = f[0]
        it = f[1]
        julian = int( TL_orig.timelist[it].strftime("%j") )
        ncIN = NC.netcdf_file(filename,'r')
        CHL_IN = ncIN.variables['CHL'].data[0,:,:]
        Lon    = ncIN.variables['Lon'].data
        Lat    = ncIN.variables['Lat'].data
        
        CHL_IN[581:,164:] = fillValue # BLACK SEA
        counter_refNAN=0
        counter_elim = 0
        
        CHL_OUT = np.zeros_like(CHL_IN)
        for i in range(jpi):
            for j in range(jpj):
                if CHL_IN[j,i] == fillValue:
                    CHL_OUT[j,i] = fillValue
                    counter_refNAN = counter_refNAN +1
                else:
                    if np.abs(CHL_IN[j,i] - DAILY_REF_MEAN[j,i]) < DAILY_REF_STD[j,i] *2.0:
                        CHL_OUT[j,i] = CHL_IN[j,i]
                    else:
                        CHL_OUT[j,i] = fillValue
                        counter_elim = counter_elim + 1
                    
        
         
    





