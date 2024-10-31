import argparse
from bitsea.utilities.argparse_types import some_among, date_from_str
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
    parser.add_argument(   '--varnames', '-v',
                                type = some_among(['CHL','KD490','DIATO','NANO','PICO', 'DINO','RRS412','RRS443','RRS490','RRS510','RRS555','RRS670']),
                                required = True
                                )

    parser.add_argument(   '--timeaverage', '-t',
                                type = str,
                                required = True,
                                choices = ['monthly','weekly_tuesday','weekly_friday','weekly_monday','weekly_thursday','tendays'],
                                help = ''' Name of the mesh of sat ORIG and used to dump checked data.'''
                                )
    parser.add_argument(   '--ignore-after', '-a',
                                type = date_from_str,
                                required = False,
                                default = None,
                                help = "Ignore all input files with dates later than the one submitted here"
                                )
    parser.add_argument(   '--force', '-f',
                                action='store_true',
                                help = """Overwrite existing variables in files
                                """)
    parser.add_argument(   '--serial',
                                action='store_true',
                                help = '''Do not use mpi''')
    return parser.parse_args()

args = argument()


from bitsea.commons.Timelist import TimeList
from bitsea.commons.time_interval import TimeInterval
from datetime import datetime
from bitsea.postproc import masks
import numpy as np
import os
import bitsea.Sat.SatManager as Sat
from bitsea.commons.utils import addsep
from bitsea.utilities.mpi_serial_interface import get_mpi_communicator

if not args.serial:
    import mpi4py.MPI


comm = get_mpi_communicator()
rank  = comm.Get_rank()
nranks =comm.size

CHECKDIR = addsep(args.checkdir)
OUTDIR   = addsep(args.outdir)
maskSat = getattr(masks,args.mesh)




Timestart = datetime.strptime("19500101", "%Y%m%d")
Time__end = datetime.strptime("20500101", "%Y%m%d")

if args.ignore_after is not None:
    Time__end = args.ignore_after

TI = TimeInterval.fromdatetimes(Timestart, Time__end)
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

    ii, w = TLCheck.select(req)
    nFiles = len(ii)
    if nFiles < 3 : 
        print(req, "less than 3 files - Skipping average generation")
        print(" ".join([TLCheck.Timelist[k].strftime("%Y%m%d") for k in ii]))
        continue

    for varname in args.varnames:
        writing_mode = Sat.writing_mode(outfile)
        if os.path.exists(outfile):
            try:
                dateweek_string = Sat.read_variable_attribute(outfile,varname,'average_of')
            except (AttributeError,IndexError):
                dateweek_string = ""
            nDates=len(dateweek_string.split(","))


            if nFiles>nDates:
                pass #print('Not skipping ' + req.string + " for " + varname)
            else:
                condition_to_write = not Sat.exist_valid_variable(varname,outfile)
                if args.force: condition_to_write=True
                if not condition_to_write: continue


        print(outfile, varname, flush=True)
        dateweek = []


        M = np.zeros((nFiles,jpj,jpi),np.float32)
        for iFrame, j in enumerate(ii):
            inputfile = TLCheck.filelist[j]
            VALUES = Sat.readfromfile(inputfile, varname)
            M[iFrame,:,:] = VALUES
            idate = TLCheck.Timelist[j]
            date8 = idate.strftime('%Y%m%d')
            dateweek.append(date8)
        if varname == 'KD490':
            OUT = Sat.averager(M)
        else:
            OUT = Sat.logAverager(M)
        dateweek_string=','.join(dateweek)
        var_attributes={'average_of':dateweek_string}
        Sat.dumpGenericfile(outfile, OUT, varname, mesh=maskSat, mode=writing_mode, var_attributes=var_attributes)

