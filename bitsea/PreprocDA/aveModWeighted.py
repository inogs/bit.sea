import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generic averager for mod surface chl
    Files with dates used for the average provided (dirdates).
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(   '--indir', '-i',
                                type = str,
                                required = True,
                                help = ''' Model directory with 3D Chl'''

                                )
    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                help = ''' OUT average dates sat directory'''

                                )


    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help = ''' maskfile'''
                                )

    parser.add_argument(   '--timeaverage', '-t',
                                type = str,
                                required = True,
                                choices = ['weekly_tuesday','weekly_friday','weekly_monday','weekly_thursday'],
                                help = ''' Type of average (weighted only on weekly)'''
                                )

    parser.add_argument(   '--stdweight', '-s',
                                type = str,
                                required = True,
                                help = ''' Std for Gaussian weight (n. days of time correlation'''
                                )

    return parser.parse_args()

args = argument()


import numpy as np
from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
from commons.mask import Mask
from commons.utils import addsep
from commons.dataextractor import DataExtractor
import Sat.SatManager as Sat
from postproc import masks
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

INDIR = addsep(args.indir)
OUTDIR   = addsep(args.outdir)
std = float(args.stdweight)
TheMask = Mask(args.maskfile)


reset = True
chlnm = 'Chla'

Timestart="19500101"
Time__end="20500101"
TI = TimeInterval(Timestart,Time__end,"%Y%m%d")
TLchl = TimeList.fromfilenames(TI, INDIR,'ave' + '*' + chlnm + '.nc', \
                               prefix='ave.')
#suffix = os.path.basename(TLchl.filelist[0])[8:]


if args.timeaverage == 'weekly_tuesday' : TIME_reqs=TLchl.getWeeklyList(2)
if args.timeaverage == 'weekly_friday'  : TIME_reqs=TLchl.getWeeklyList(5)
if args.timeaverage == 'weekly_monday'  : TIME_reqs=TLchl.getWeeklyList(1)
if args.timeaverage == 'weekly_thursday': TIME_reqs=TLchl.getWeeklyList(4)

_,jpj,jpi = TheMask.shape
masknan = TheMask.mask[0,:,:]==False

if (jpj==380) & (jpi==1085):
   maskSat = getattr(masks,'Mesh24')
else:
   print 'mask sat to be added'
   import sys
   sys.exit(0)

counter = 0
MySize = len(TIME_reqs[rank::nranks])

fillValue = -999

for req in TIME_reqs[rank::nranks]:
    counter = counter + 1

    outfile = 'ave.' + req.string + '-00:00:00.' + chlnm + '.nc'
    outpathfile = OUTDIR + outfile

    ii, w = TLchl.selectWeeklyGaussWeights(req,std)
    nFiles = len(ii)

    print outfile
    dateweek = []
    if nFiles < 3 : 
        print req
        print "less than 3 files"
        continue
    M = np.zeros((nFiles,jpj,jpi),np.float32)
    for iFrame, j in enumerate(ii):
        inputfile = TLchl.filelist[j]
        De = DataExtractor(TheMask,filename=inputfile,varname=chlnm,dimvar=3)
        CHL = De.filled_values[0,:,:].copy()
        CHL[masknan] = fillValue
        M[iFrame,:,:] = CHL
    CHL_OUT = Sat.WeightedLogAverager(M,w)
    Sat.dumpGenericNativefile(outpathfile, CHL_OUT, varname=chlnm, mesh=maskSat)




    print "\trequest ", counter, " of ", MySize, " done by rank ", rank
