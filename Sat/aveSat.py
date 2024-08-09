import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generic averager for sat files.
    It works with one mesh, without performing interpolations.
    Empty dirdates is provided when the number of input daily files is lesser than 3.
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

    parser.add_argument(   '--mesh', '-m',
                                type = str,
                                required = True,
                                choices = ['SatOrigMesh','V4mesh','V1mesh','KD490mesh','SAT1km_mesh', 'Mesh24', 'Mesh4'],
                                help = ''' Name of the mesh of sat ORIG and used to dump checked data.'''
                                )
    parser.add_argument(   '--varname', '-v',
                                type = str,
                                required = True,
                                choices = ['CHL','KD490','DIATO','NANO','PICO', 'DINO','RRS412','RRS443','RRS490','RRS510','RRS555','RRS670']
                                )

    parser.add_argument(   '--timeaverage', '-t',
                                type = str,
                                required = True,
                                choices = ['monthly','weekly_tuesday','weekly_friday','weekly_monday','weekly_thursday','tendays'],
                                help = ''' Name of the mesh of sat ORIG and used to dump checked data.'''
                                )
    parser.add_argument(   '--force', '-f',
                                action='store_true',
                                help = """Overwrite existing variables in files
                                """)
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
maskSat = getattr(masks,args.mesh)



Timestart="19500101"
Time__end="20500101"
TI = TimeInterval(Timestart,Time__end,"%Y%m%d")
TLCheck = TimeList.fromfilenames(TI, CHECKDIR,"*.nc",prefix='',dateformat='%Y%m%d')
suffix = os.path.basename(TLCheck.filelist[0])[8:]


if args.timeaverage == 'monthly'        : TIME_reqs=TLCheck.getMonthlist()
if args.timeaverage == 'weekly_tuesday' : TIME_reqs=TLCheck.getWeeklyList(2)
if args.timeaverage == 'weekly_friday'  : TIME_reqs=TLCheck.getWeeklyList(5)
if args.timeaverage == 'weekly_monday'  : TIME_reqs=TLCheck.getWeeklyList(1)
if args.timeaverage == 'weekly_thursday': TIME_reqs=TLCheck.getWeeklyList(4)
if args.timeaverage == 'tendays'        : TIME_reqs=TLCheck.getSpecificIntervalList(10,"19971001-12:00:00")

jpi = maskSat.jpi
jpj = maskSat.jpj



for req in TIME_reqs[rank::nranks]:

    outfile = OUTDIR + req.string + suffix
    writing_mode = Sat.writing_mode(outfile)

    ii, w = TLCheck.select(req)
    nFiles = len(ii)

    if os.path.exists(outfile):
        try:
            dateweek_string = Sat.read_variable_attribute(outfile,args.varname,'average_of')
        except (AttributeError,IndexError):
            dateweek_string = ""
        nDates=len(dateweek_string.split(","))

    
        if nFiles>nDates:
            print('Not skipping ' + req.string)
        else:
            condition_to_write = not Sat.exist_valid_variable(args.varname,outfile)
            if args.force: condition_to_write=True
            if not condition_to_write: continue


    print(outfile, flush=True)
    dateweek = []
    if nFiles < 3 : 
        print(req, "less than 3 files - Skipping average generation")
        print(" ".join([TLCheck.Timelist[k].strftime("%Y%m%d") for k in ii]))
        continue

    M = np.zeros((nFiles,jpj,jpi),np.float32)
    for iFrame, j in enumerate(ii):
        inputfile = TLCheck.filelist[j]
        VALUES = Sat.readfromfile(inputfile, args.varname)
        M[iFrame,:,:] = VALUES
        idate = TLCheck.Timelist[j]
        date8 = idate.strftime('%Y%m%d')
        dateweek.append(date8)
    if args.varname == 'KD490':
        OUT = Sat.averager(M)
    else:
        OUT = Sat.logAverager(M)
    dateweek_string=','.join(dateweek)
    var_attributes={'average_of':dateweek_string}
    Sat.dumpGenericfile(outfile, OUT, args.varname, mesh=maskSat, mode=writing_mode, var_attributes=var_attributes)

