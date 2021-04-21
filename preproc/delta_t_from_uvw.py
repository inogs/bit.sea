import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates a txt file for each forcing time
    Example of name:
       DeltaT_19990101-00:00:00.txt
    Content:
      900  703.634   37 138 298
    ogstm.xx reads the first number, an integer which can be 600, 450, 360, 300
    ''')
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help = '/gpfs/work/OGS_prod_0/OPA/V5C/devel/wrkdir/2/MODEL/FORCINGS/')
    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                help = '/some/path/')
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help = '/gpfs/scratch/userexternal/gbolzon0/OPEN_BOUNDARY/FORCINGS_CHECK/meshmask_INGV.nc')

    return parser.parse_args()

args = argument()

from commons.mask import Mask
from commons.dataextractor import DataExtractor
import numpy as np
from commons.Timelist import TimeList
from commons.utils import addsep

try:
    from mpi4py import MPI
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    nranks =comm.size
except:
    rank   = 0
    nranks = 1

INPUTDIR=addsep(args.inputdir)
OUTPUTDIR=addsep(args.outdir)
TheMask=Mask(args.maskfile)

jpk, jpj, jpi = TheMask.shape
E1T=np.zeros((jpk,jpj,jpi),np.float32)
E2T=np.zeros((jpk,jpj,jpi),np.float32)
for k in range(jpk):
    E1T[k,:,:] = TheMask.e1t
for k in range(jpk):
    E2T[k,:,:] = TheMask.e2t


def impose_deltat(deltaT):
    if deltaT > 575: return 600
    if deltaT > 435: return 450
    if deltaT > 345: return 360
    return 300



TL=TimeList.fromfilenames(None, INPUTDIR, "U*nc", prefix="U",hour=0)
eps=1.e-08
Cmax=1.0




nFrames=TL.nTimes
mydtype=[('Imposed_deltaT',np.int), ('deltaT',np.float32),('K',np.int),('J',np.int),('I',np.int)]

FRAMES=range(nFrames)

for iframe in FRAMES[rank::nranks]:
    DELTAT=np.zeros((1,),dtype=mydtype)
    timestr = TL.Timelist[iframe].strftime("%Y%m%d-%H:%M:%S")
    print timestr
    filenameU=INPUTDIR + "U" + timestr + ".nc"
    filenameV=INPUTDIR + "V" + timestr + ".nc"
    filenameW=INPUTDIR + "W" + timestr + ".nc"

    U=np.abs(DataExtractor(TheMask,filenameU,"vozocrtx").values)
    V=np.abs(DataExtractor(TheMask,filenameV,"vomecrty").values)
    W=np.abs(DataExtractor(TheMask,filenameW,"vovecrtz").values)
    
    U[U==0]=eps
    V[V==0]=eps
    W[W==0]=eps
    U[U>1.e+19]=eps
    V[V>1.e+19]=eps
    W[W>1.e+19]=eps
    Fact = U/E1T + V/E2T + W/TheMask.e3t
    deltat = Cmax/Fact
    K,J,I = np.nonzero(deltat==deltat.min())
    calculated_deltat = deltat.min()
    DELTAT[0]['Imposed_deltaT']= impose_deltat(calculated_deltat)
    DELTAT[0]['deltaT']=calculated_deltat
    DELTAT[0]['K']=K[0]
    DELTAT[0]['J']=J[0]
    DELTAT[0]['I']=I[0]
    outfile=OUTPUTDIR + "DeltaT_" + timestr + ".txt"
    np.savetxt(outfile, DELTAT,fmt="%5d %10.3f %d %d %d")
    


