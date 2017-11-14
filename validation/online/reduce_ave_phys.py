import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Creates T,S files archive, reading from ave phys.
    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(   '--inputdir','-i',
                                type = str,
                                required = True,
                                help = 'input dir validation tmp')
    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                default = None,
                                required = True,
                                help = "profile dir")
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                default = None,
                                required = True,
                                help = ''' Path of maskfile''')

    return parser.parse_args()

args = argument()
import numpy as np
from commons.utils import addsep
from commons.dataextractor import DataExtractor
from commons.mask import Mask
from commons.Timelist import TimeInterval, TimeList
from commons import netcdf4
try:
    from mpi4py import MPI
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    nranks =comm.size
except:
    rank   = 0
    nranks = 1


INPUTDIR  = addsep(args.inputdir)
OUTPUTDIR = addsep(args.outdir)
TheMask    = Mask(args.maskfile)

TI = TimeInterval("1950","2050",'%Y')
TL = TimeList.fromfilenames(TI, INPUTDIR, "ave*phys.nc")
nTimes = TL.nTimes
VARLIST=['vosaline','votemper']
nvars = len(VARLIST)

PROCESSES = np.arange(nTimes * nvars)
for ip in PROCESSES[rank::nranks]:
    (itime, ivar) = divmod(ip,nvars)
    var = VARLIST[ivar]
    
    phys_file = INPUTDIR + TL.Timelist[itime].strftime("ave.%Y%m%d-12:00:00.phys.nc")
    outfile  = OUTPUTDIR + TL.Timelist[itime].strftime("ave.%Y%m%d-12:00:00.") + var + ".nc"
    print "rank = ", rank, var, phys_file, outfile
    T = DataExtractor(TheMask,phys_file,var).values
    netcdf4.write_3d_file(T, var, outfile, TheMask)
