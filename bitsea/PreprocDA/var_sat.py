import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Computes the Monthly surface variance of the dayF-days averages
    The inputs file are within INDIR directory and the outputs will be placed
    within OUTDIR

    The mesh must be consistent with the data contained in INDIR
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(   '--outdir', '-o',
                            type = str,
                            required =True,
                            default = "./",
                            help = ''' Outputdir'''
                            )

    parser.add_argument(   '--indir', '-i',
                            type = str,
                            required = True,
                            help = 'Input map sat'
                            )

    parser.add_argument(   '--mesh', '-m',
                            type = str,
                            required = True,
                            help = 'Mesh of sat map resolution'
                            )

    parser.add_argument(   '--limstd', '-s',
                            type = str,
                            required = True,
                            help = 'std used for outlyers exclusion (used only for output files names)'
                            )

    parser.add_argument(   '--dayf', '-f',
                            type = str,
                            required = True,
                            help = 'number of days sat map frequency (used only for output files names)'
                            )

    return parser.parse_args()

args = argument()
import numpy as np
from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
from commons.timerequestors import Clim_month
from commons.utils import addsep
from postproc import masks
import Sat.SatManager as Sat
import netCDF4


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



limstd = args.limstd
dayF = args.dayf
MyMesh = getattr(masks,args.mesh)
INDIR = addsep(args.indir)
OUTDIR = addsep(args.outdir)



Timestart="19500101"
Time__end="20500101"

TI = TimeInterval(Timestart,Time__end,"%Y%m%d")
TLCheck = TimeList.fromfilenames(TI, INDIR,"*.nc",prefix='',dateformat='%Y%m%d')
MONTHLY_reqs=[Clim_month(i) for i in range(1,13)]

print "number of requests: ", len(MONTHLY_reqs), " inputFreq: ", TLCheck.inputFrequency

jpi = MyMesh.jpi
jpj = MyMesh.jpj

Chl = np.zeros((12,jpj,jpi), dtype=float)
ChlSquare = np.zeros((12,jpj,jpi), dtype=float)

# get fillValue assuming that all files in the
# folder have the same fillValue
ncIN = netCDF4.Dataset(TLCheck.filelist[0],'r')
varObj = ncIN.variables['CHL']
fillValue = varObj.getncattr('_FillValue')
ncIN.close()

for req in MONTHLY_reqs[rank::nranks]:
    ii, w = TLCheck.select(req)

    nFiles = len(ii)
    M  = np.zeros((nFiles,jpj,jpi),np.float32)
    M2 = np.zeros((nFiles,jpj,jpi),np.float32)
    
    for iFrame, j in enumerate(ii):
        inputfile = TLCheck.filelist[j]
        CHL = Sat.readfromfile(inputfile,'CHL')
        TmpMat = np.zeros(CHL.shape)
        M[iFrame,:,:]  = CHL

        MyMask = np.where(CHL==fillValue)
        TmpMat = CHL*CHL
        TmpMat[MyMask] = fillValue
        M2[iFrame,:,:] = TmpMat
    
    MonthIndex = req.month-1
    Chl[MonthIndex,:,:] = Sat.WeightedAverager(M,w)
    ChlSquare[MonthIndex,:,:] = Sat.WeightedAverager(M2,w)

#for ii in range(0,12):
    # computing variance
    Var2D = ChlSquare[MonthIndex,:,:] - Chl[MonthIndex,:,:]*Chl[MonthIndex,:,:]
    
    # filling missing points
    MyMask = np.where((ChlSquare[MonthIndex,:,:] == fillValue) | (Var2D<0))
    Var2D[MyMask]  = fillValue
    
    # saving results
    MonthStr = '%02d'%(MonthIndex+1)
    if OUTDIR[-1]=="/":
      fname = 'var2Dsat.CCI.F'
    else:
      fname = '/var2Dsat.CCI.F'
    
    fname += dayF+'.'+limstd+'.'+MonthStr+".nc"
    filename = OUTDIR+fname
    print "\tsaving ", fname
    Sat.dumpGenericNativefile(filename,Var2D,varname='variance',mesh=MyMesh)

print "\n\tDone :)\n"

