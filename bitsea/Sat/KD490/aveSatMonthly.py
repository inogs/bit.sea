from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
import numpy as np
import os
import Sat.SatManager as Sat

CHECKDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/KD490/DAILY/CHECKED/"
MONTHLYDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/KD490/MONTHLY/ORIGMESH/"

reset = False

Timestart="19990101"
Time__end="20160101"
TI = TimeInterval(Timestart,Time__end,"%Y%m%d")
TLCheck = TimeList.fromfilenames(TI, CHECKDIR,"*.nc",prefix='',dateformat='%Y%m%d')
#IonamesFile = '../../postproc/IOnames_sat.xml'
#IOname = IOnames.IOnames(IonamesFile)
suffix = os.path.basename(TLCheck.filelist[0])[8:]
MONTHLY_reqs=TLCheck.getMonthlist()

jpi = Sat.masks.KD490mesh.jpi
jpj = Sat.masks.KD490mesh.jpj


for req in MONTHLY_reqs:
    outfile = req.string + suffix
    outpathfile = MONTHLYDIR + outfile
    conditionToSkip = (os.path.exists(outpathfile)) and (not reset)

    if conditionToSkip: continue

    print outfile
    ii, w = TLCheck.select(req)
    nFiles = len(ii)
    M = np.zeros((nFiles,jpj,jpi),np.float32)
    for iFrame, j in enumerate(ii):
        inputfile = TLCheck.filelist[j]
        CHL = Sat.readfromfile(inputfile,'KD490')
        M[iFrame,:,:] = CHL
    Kext_OUT = Sat.averager(M)
    Sat.dumpGenericNativefile(outpathfile, Kext_OUT)
