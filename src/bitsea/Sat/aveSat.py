import argparse
from bitsea.utilities.argparse_types import some_among, date_from_str, existing_dir_path

def argument():
    parser = argparse.ArgumentParser(description = '''
    Generic averager for sat files.
    It works with one mesh, without performing interpolations.
    Empty dirdates is provided when the number of input daily files is lesser than 3.
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(   '--checkdir', '-i',
                                type = existing_dir_path,
                                required = True,
                                help = ''' CHECKED sat directory, e.g. /gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/ONLINE/SAT/MODIS/DAILY/CHECKED/'''

                                )
    parser.add_argument(   '--outdir', '-o',
                                type = existing_dir_path,
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

from pathlib import Path
from typing import List, Iterable
from bitsea.commons.Timelist import TimeList
from bitsea.commons.time_interval import TimeInterval
from datetime import datetime
from bitsea.postproc import masks
import numpy as np
import bitsea.Sat.SatManager as Sat
from sys import exit

from bitsea.utilities.mpi_serial_interface import get_mpi_communicator


def aveSat(*,
inputfiles: Iterable[Path],
datetimes: List[str],
outputdir: Path,
mesh: str,
varnames = List[str],
serial : bool = False,
outfrequency : str = "weekly",
week_day : int = 1,
force : bool,
)-> List[Path]:

    if not serial:
        import mpi4py.MPI


    comm = get_mpi_communicator()
    rank  = comm.Get_rank()
    nranks =comm.size

    OUTDIR   = outputdir
    maskSat = getattr(masks,mesh)


    suffix = inputfiles[0].name[8:]
    TLCheck = TimeList(datetimes)
    if outfrequency == 'monthly'    : TIME_reqs=TLCheck.getMonthlist()
    if outfrequency == "weekly"     : TIME_reqs=TLCheck.getWeeklyList(week_day)
    if outfrequency == "ten_days"   : TIME_reqs=TLCheck.getSpecificIntervalList(10,"19971001-12:00:00")


    jpi = maskSat.jpi
    jpj = maskSat.jpj
    averaged_filelist = [ OUTDIR / (req.string + suffix) for req in TIME_reqs ]



    for req in TIME_reqs[rank::nranks]:

        outfile = OUTDIR / (req.string + suffix)

        ii, w = TLCheck.select(req)
        nFiles = len(ii)
        if nFiles < 3 :
            print(req, "less than 3 files - Skipping average generation")
            print(" ".join([TLCheck.Timelist[k].strftime("%Y%m%d") for k in ii]))
            continue

        for varname in varnames:
            writing_mode = Sat.writing_mode(outfile)
            if outfile.exists():
                try:
                    dateweek_string = Sat.read_variable_attribute(outfile,varname,'average_of')
                except (AttributeError,IndexError):
                    dateweek_string = ""
                nDates=len(dateweek_string.split(","))


                if nFiles>nDates:
                    pass #print('Not skipping ' + req.string + " for " + varname)
                else:
                    condition_to_write = not Sat.exist_valid_variable(varname,outfile)
                    if force: condition_to_write=True
                    if not condition_to_write: continue


            print(outfile, varname, flush=True)
            dateweek = []


            M = np.zeros((nFiles,jpj,jpi),np.float32)
            for iFrame, j in enumerate(ii):
                inputfile = inputfiles[j]
                VALUES = Sat.readfromfile(inputfile, varname)
                M[iFrame,:,:] = VALUES
                idate = datetimes[j]
                date8 = idate.strftime('%Y%m%d')
                dateweek.append(date8)
            if varname == 'KD490':
                OUT = Sat.averager(M)
            else:
                OUT = Sat.logAverager(M)
            dateweek_string=','.join(dateweek)
            var_attributes={'average_of':dateweek_string}
            Sat.dumpGenericfile(outfile, OUT, varname, mesh=maskSat, mode=writing_mode, var_attributes=var_attributes)
    return averaged_filelist


def main():
    args = argument()

    Timestart = datetime.strptime("19500101", "%Y%m%d")
    Time__end = datetime.strptime("20500101", "%Y%m%d")

    if args.ignore_after is not None:
        Time__end = args.ignore_after

    TI = TimeInterval.fromdatetimes(Timestart, Time__end)
    TL = TimeList.fromfilenames(TI, args.checkdir,"*.nc",prefix='',dateformat='%Y%m%d')

    outfrequency = 'weekly'
    weekday = None
    if args.timeaverage == 'monthly'        : outfrequency = 'monthly'
    if args.timeaverage == 'weekly_tuesday' : weekday=2
    if args.timeaverage == 'weekly_friday'  : weekday=5
    if args.timeaverage == 'weekly_monday'  : weekday=1
    if args.timeaverage == 'weekly_thursday': weekday=4
    if args.timeaverage == 'tendays'        : outfrequency = 'tendays'

    aveSat(
        inputfiles=TL.filelist,
        datetimes = TL.Timelist,
        outputdir=args.outdir,
        mesh=args.mesh,
        varnames=args.varnames,
        outfrequency=outfrequency,
        week_day = weekday,
        force = args.force,
        serial=args.serial,
    )
    return 0


if __name__ == "__main__":
    exit(main())