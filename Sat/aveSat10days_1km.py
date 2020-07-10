from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
from commons import IOnames
import numpy as np
import os
import Sat.SatManager as Sat
try:
    from mpi4py import MPI
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    nranks =comm.size
    isParallel = True
except:
    rank   = 0
    nranks = 1
    isParallel = False

CHECKDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/CCI_1km/DAILY/CHECKED/"
OUTDIR="/pico/home/userexternal/pdicerbo/WorkDir/AveSat24/Checked_10Days_Sat1km/"

reset = False

Timestart="19500101"
Time__end="20500101"
TI = TimeInterval(Timestart,Time__end,"%Y%m%d")
TLCheck = TimeList.fromfilenames(TI, CHECKDIR,"*.nc",prefix='',dateformat='%Y%m%d')
IonamesFile = '../../postproc/IOnames_sat.xml'
IOname = IOnames.IOnames(IonamesFile)

TenDaysReqs = TLCheck.getSpecificIntervalList(deltastr='days=10', starttime="19500101-12:00:00")

jpi = Sat.OneKmMesh.jpi
jpj = Sat.OneKmMesh.jpj

counter = 0
MySize = len(TenDaysReqs[rank::nranks])

for req in TenDaysReqs[rank::nranks]:
    counter = counter + 1
    
    outfile = req.string + IOname.Output.suffix + ".nc"
    outpathfile = OUTDIR + outfile
    conditionToSkip = (os.path.exists(outpathfile)) and (not reset)

    if conditionToSkip: continue

    print outfile
    ii, w = TLCheck.select(req)
    nFiles = len(ii)
    if nFiles < 3 : 
        print req
        print "less than 3 files"
        continue
    M = np.zeros((nFiles,jpj,jpi),np.float32)
    for iFrame, j in enumerate(ii):
        inputfile = TLCheck.filelist[j]
        CHL = Sat.readfromfile(inputfile)
        M[iFrame,:,:] = CHL
    CHL_OUT = Sat.logAverager(M)
    Sat.dumpGenericNativefile(outpathfile, CHL_OUT, varname='CHL', mesh=Sat.OneKmMesh)

    print "\trequest ", counter, " of ", MySize, " done by rank ", rank