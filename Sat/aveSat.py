import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generic averager for sat files.
    It works with one mesh, without performing interpolations.
    Files with dates used for the average provided (dirdates).
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(   '--checkdir', '-i',
                                type = str,
                                required = True,
                                help = ''' CHECKED sat directory, e.g. /gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MODIS/DAILY/CHECKED/'''

                                )
    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                help = ''' OUT average dates sat directory'''

                                )

    parser.add_argument(   '--dirdates', '-d',
                                type = str,
                                required = True,
                                help = ''' OUT dates sat directory'''

                                )

    parser.add_argument(   '--mesh', '-m',
                                type = str,
                                required = True,
                                choices = ['SatOrigMesh','V4mesh','V1mesh','KD490mesh','SAT1km_mesh', 'Mesh24'],
                                help = ''' Name of the mesh of sat ORIG and used to dump checked data.'''
                                )

    parser.add_argument(   '--timeaverage', '-t',
                                type = str,
                                required = True,
                                choices = ['monthly','weekly_tuesday','weekly_friday'],
                                help = ''' Name of the mesh of sat ORIG and used to dump checked data.'''
                                )
    return parser.parse_args()

args = argument()


from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
from postproc import masks
import numpy as np
import os
import Sat.SatManager as Sat
from commons.utils import addsep
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

CHECKDIR = addsep(args.checkdir)
OUTDIR   = addsep(args.outdir)
DIRDATES = addsep(args.dirdates)
maskSat = getattr(masks,args.mesh)

reset = True

Timestart="19500101"
Time__end="20500101"
TI = TimeInterval(Timestart,Time__end,"%Y%m%d")
TLCheck = TimeList.fromfilenames(TI, CHECKDIR,"*.nc",prefix='',dateformat='%Y%m%d')
suffix = os.path.basename(TLCheck.filelist[0])[8:]

if args.timeaverage == 'monthly'        : TIME_reqs=TLCheck.getMonthlist()
if args.timeaverage == 'weekly_tuesday' : TIME_reqs=TLCheck.getWeeklyList(2)
if args.timeaverage == 'weekly_friday'  : TIME_reqs=TLCheck.getWeeklyList(5)

jpi = maskSat.jpi
jpj = maskSat.jpj

counter = 0
MySize = len(TIME_reqs[rank::nranks])

for req in TIME_reqs[rank::nranks]:
    counter = counter + 1

    outfile = req.string + suffix
    outpathfile = OUTDIR + outfile
    conditionToSkip = (os.path.exists(outpathfile)) and (not reset)

    if conditionToSkip: continue

    print outfile
    ii, w = TLCheck.select(req)
    nFiles = len(ii)
    dateweek = []
    if nFiles < 3 : 
        print req
        print "less than 3 files"
        filedates = DIRDATES + req.string + 'weekdates.txt'
        print(filedates)
        np.savetxt(filedates,dateweek,fmt='%s')
        continue
    M = np.zeros((nFiles,jpj,jpi),np.float32)
    for iFrame, j in enumerate(ii):
        inputfile = TLCheck.filelist[j]
        CHL = Sat.readfromfile(inputfile)
        M[iFrame,:,:] = CHL
        idate = TLCheck.Timelist[j]
        date8 = idate.strftime('%Y%m%d')
        dateweek.append(date8)
    CHL_OUT = Sat.logAverager(M)
    Sat.dumpGenericNativefile(outpathfile, CHL_OUT, varname='CHL', mesh=maskSat)


    filedates = DIRDATES + req.string + 'weekdates.txt'
    print(filedates)
    np.savetxt(filedates,dateweek,fmt='%s')

    print "\trequest ", counter, " of ", MySize, " done by rank ", rank
