import netCDF4
import numpy as np
from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
import dom_dec
import time_manager
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

maskfile="CHL_1km_meshmask.nc"
INPUTDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/CCI_1km/DAILY/ORIG"

D=netCDF4.Dataset(maskfile,'r')
tmask_glo = np.array(D.variables['tmask']).astype(np.bool)
D.close()

jpjglo,jpiglo = tmask_glo.shape

nproc_i = 4
nproc_j = 4
PROCESSES = np.arange(nproc_i*nproc_j)
valid_max=15.0


for ip in PROCESSES[rank::nranks]:
    (jproc, iproc) = divmod(ip,nproc_i)
    inputfile = "CHL_raw_data_%d_%d.npy" %(jproc,iproc)
    outputfile= "CHL_clim_%d_%d" %(jproc,iproc)
    I_start,I_end, J_start, J_end = dom_dec.dom_dec(iproc, jproc, jpiglo, jpjglo, nproc_i, nproc_j)
    tmask = tmask_glo[J_start:J_end, I_start:I_end]
    print "rank ", rank, "processes " , ip, outputfile, iproc, jproc, I_start,I_end, J_start, J_end
    Nwaterpoints = tmask.sum()

    print inputfile
    print "Nwaterpoints", Nwaterpoints, ip
    if (Nwaterpoints ==0) : continue


    RAW_DATA = np.load(inputfile)



    WaterPoints = np.zeros_like(tmask, dtype=np.int32)
    WaterPoints[tmask] = np.arange(Nwaterpoints)
    print "Nwaterpoints", Nwaterpoints
    jpi = I_end-I_start
    jpj = J_end-J_start
    CLIM = np.zeros((365, jpj,jpi), dtype=[('NUMB',np.int32), ('MEAN',np.float32),('STD',np.float32)])
    for julian in range(365):
        print "rank ", rank, " julian day " ,julian
        II, filelist=time_manager.getfilelist(julian)
        raw_data_julian=RAW_DATA[:,II]
        for ji in range(jpi):
            for jj in range(jpj):
                if tmask[jj,ji]:
                    jw = WaterPoints[jj,ji];
                    PILAloc = raw_data_julian[jw,:];
                    valid_points = (PILAloc > -999) & (PILAloc < valid_max)
                    if valid_points.sum() < 5: continue
                    pilarray= PILAloc[valid_points];

                    m = pilarray.mean()
                    s = pilarray.std()
                    outsiders = pilarray > m+3.0*s
                    pilarray=pilarray[~outsiders]
                    n = len(pilarray)
                    CLIM['NUMB'][julian, jj,ji] = n
                    if n>4:
                        m = pilarray.mean()
                        CLIM['MEAN'][julian,jj,ji] = m
                        CLIM['STD' ][julian,jj,ji] =np.sqrt( ((pilarray-m)**2).sum()/(n-1)  )

    np.save(outputfile,CLIM)
        
