from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
from commons import IOnames
import numpy as np
import os
import Sat.SatManager as Sat

CHECKDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/CCI/DAILY/CHECKED/"
MONTHLYDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/CCI/MONTHLY_V4/"

reset = False

Timestart="19990102"
Time__end="20160101"
TI = TimeInterval(Timestart,Time__end,"%Y%m%d")
TLCheck = TimeList.fromfilenames(TI, CHECKDIR,"*.nc",prefix='',dateformat='%Y%m%d')
IonamesFile = '../../postproc/IOnames_sat.xml'
IOname = IOnames.IOnames(IonamesFile)

MONTHLY_reqs=TLCheck.getMonthlist()

jpi = Sat.NativeMesh.jpi
jpj = Sat.NativeMesh.jpj


for req in MONTHLY_reqs:
    outfile = req.string + IOname.Output.suffix + ".nc"
    outpathfile = MONTHLYDIR + outfile
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
    Sat.dumpV4file(outpathfile, CHL_OUT)

