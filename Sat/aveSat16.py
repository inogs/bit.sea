import postproc.Timelist as Timelist
import postproc.IOnames as IOnames
import numpy as np
import os
import SatManager as Sat

CHECKDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MODIS/DAILY/CHECKED/"
WEEKLYDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MODIS/WEEKLY/"

reset = False

Timestart="19500101"
Time__end="20500101"
IonamesFile = '../postproc/IOnames_sat.xml'
TLCheck = Timelist.TimeList(Timestart,Time__end, CHECKDIR,"*.nc",IonamesFile)
IOname = IOnames.IOnames(IonamesFile)

WEEK_reqs=TLCheck.getWeeklyList(2)

jpi = Sat.NativeMesh.jpi
jpj = Sat.NativeMesh.jpj


for req in WEEK_reqs:
    outfile = req.string + IOname.Output.suffix + ".nc"
    outpathfile = WEEKLYDIR + outfile
    conditionToSkip = (os.path.exists(outpathfile)) and (not reset)

    if conditionToSkip: continue

    print outfile
    ii, w = TLCheck.select(req)
    nFiles = len(ii)
    M = np.zeros((nFiles,jpj,jpi),np.float32)
    for iFrame, j in enumerate(ii):
        inputfile = TLCheck.filelist[j]
        CHL = Sat.readfromfile(inputfile)
        M[iFrame,:,:] = CHL
    CHL_OUT = Sat.logAverager(M)
    Sat.dumpfile(outpathfile, CHL_OUT)

