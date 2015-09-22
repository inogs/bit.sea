import postproc
from postproc import Timelist
from postproc import IOnames
import numpy as np

import SatManager as Sat

CHECKDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MODIS/DAILY/CHECKED/"
WEEKLYDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MODIS/WEEKLY/"
TLCheck = Timelist.TimeList("19500101","20500101", CHECKDIR,"*.nc",'postproc/IOnames_sat.xml')
IOname = IOnames.IOnames('postproc/IOnames_sat.xml')

WEEK_reqs=TLCheck.getWeeklyList(2)

jpi = Sat.NativeMesh.jpi
jpj = Sat.NativeMesh.jpj

for req in WEEK_reqs:
    outfile = req.string + IOname.Output.suffix + ".nc"
    ii, w = TLCheck.select(req)
    nFiles = len(ii)
    M = np.zeros((nFiles,jpj,jpi),np.float32)
    for iFrame, j in enumerate(ii):
        inputfile = TLCheck.filelist[j]
        CHL = Sat.readfromfile(inputfile)
        M[iFrame,:,:] = CHL
    
    CHL_OUT = Sat.logAverager(M)
    Sat.dumpfile(outfile, CHL_OUT)
            

