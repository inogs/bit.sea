import postproc
from postproc import Timelist
from postproc import IOnames
import glob,os
import numpy as np
import scipy.io.netcdf as NC

from sat_check import satread,satwrite, fillValue

CHECKDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MODIS/DAILY/CHECKED/"
WEEKLYDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MODIS/WEEKLY/"
TLCheck = Timelist.TimeList("19500101","20500101", CHECKDIR,"*.nc",'postproc/IOnames_sat.xml')
IOname = IOnames.IOnames('postproc/IOnames_sat.xml')

WEEK_reqs=TLCheck.getWeeklyList(2)


jpi = 733
jpj = 253

for req in WEEK_reqs:
    outfile = req.string + IOname.Output.suffix + ".nc"
    ii, w = TLCheck.select(req)
    nFiles = len(ii)
    M = np.zeros((nFiles,jpj,jpi),np.float32)
    for iFrame, j in enumerate(ii):
        inputfile = TLCheck.filelist[j]
        CHL,Lat,Lon = satread(inputfile)
        M[iFrame,:,:] = CHL
    
    CHL_OUT = np.ones((jpj,jpi),np.float32)*fillValue
    
    for i in range(jpi):
        for j in range(jpj):
            l = M[:,j,i]
            goodValues = l != fillValue
            if np.any(goodValues):
                count = goodValues.sum()
                LOGmean = np.log(l[goodValues]).sum() / count
                CHL_OUT[j,i] = np.exp(LOGmean)
    satwrite(outfile, CHL_OUT, Lon, Lat)
            

