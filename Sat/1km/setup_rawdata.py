import netCDF4
import numpy as np
from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
import Sat.SatManager as Sat
import dom_dec

try:
    from mpi4py import MPI
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    nranks =comm.size
    isParallel = True
    print "MyId = ", rank
except:
    rank   = 0
    nranks = 1
    isParallel = False

# maskfile="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/KD490/KD490_1km_meshmask.nc"
maskfile="CHL_1km_meshmask.nc"
INPUTDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/CCI_1km/DAILY/ORIG"
TI = TimeInterval("19500101","20500101","%Y%m%d")
TL = TimeList.fromfilenames(TI, INPUTDIR,"*.nc",prefix='',dateformat='%Y%m%d')
nFrames = TL.nTimes

D=netCDF4.Dataset(maskfile,'r')
tmask_glo = np.array(D.variables['tmask']).astype(np.bool)
D.close()

jpjglo,jpiglo = tmask_glo.shape

nproc_i = 4
nproc_j = 4
PROCESSES = np.arange(nproc_i*nproc_j)


for ip in PROCESSES[rank::nranks]:
    (jproc, iproc) = divmod(ip,nproc_i)
    outfile = "CHL_raw_data_%d_%d" %(jproc,iproc)

    I_start,I_end, J_start, J_end = dom_dec.dom_dec(iproc, jproc, jpiglo, jpjglo, nproc_i, nproc_j)
    tmask = tmask_glo[J_start:J_end, I_start:I_end]
    print "rank ", ip, "processes " , outfile, iproc, jproc, I_start,I_end, J_start, J_end

    Nwaterpoints = tmask.sum()
    print "Nwaterpoints", Nwaterpoints, ip
    if (Nwaterpoints ==0) : continue

    RAW_DATA = np.zeros((Nwaterpoints,nFrames),np.float32)
    for ifile, filename in enumerate(TL.filelist):
        print "reading", ip, ifile
        A=Sat.readfromfile(filename,'CHL')
        A=A[J_start:J_end, I_start:I_end]
        RAW_DATA[:,ifile] = A[tmask]

    np.save(outfile,RAW_DATA)







