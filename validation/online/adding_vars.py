import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''Adds 3D  variables to big ave files for archive and/or big files with aggregated var for 
    Avescan.py. Works in parallel.''')
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help = '/some/path/MODEL/AVE_FREQ_1/')

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                help = '''Directory containg files that will be appended
                                ''')
    parser.add_argument(   '--varname', '-v',
                                type = str,
                                required=True,
                                help = '''Name of the variable in input files, will be copied in output files as is.
                                ''')   
    parser.add_argument(   '--input_glob_pattern',"-ip",
                                type = str,
                                required=True,
                                help = 'ave*.N1p.nc, they configure the date list')
    parser.add_argument(   '--output_glob_pattern',"-op",
                                type = str,
                                required = True,
                                help = 'ave*.N1p.nc, they configure the date list')
    parser.add_argument(   '--input_prefix',
                                type = str,
                                default = "ave.",
                                help = 'The part of filename preceeding the date string')
    parser.add_argument(   '--output_prefix',
                                type = str,
                                default = "ave.",
                                help = 'The part of filename preceeding the date string')
    parser.add_argument(   '--starttime','-st',
                                type = str,
                                required = True,
                                help = 'start date in yyyymmdd format')
    parser.add_argument(   '--endtime','-et',
                                type = str,
                                required = True,
                                help = 'end date in yyyymmdd format')
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                default = None,
                                required = True,
                                help = ''' Path of maskfile''')     

    return parser.parse_args()

args = argument()

from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
from commons import netcdf3
from commons.mask import Mask
from commons.utils import addsep
from commons.dataextractor import DataExtractor
try :
    from mpi4py import MPI
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    nranks = comm.size 
except:
    rank   = 0
    nranks = 1



INPUTDIR=addsep(args.inputdir)
OUTDIR  =addsep(args.outdir)
varname = args.varname
starttime=args.starttime
end__time=args.endtime
TheMask = Mask(args.maskfile)


TI = TimeInterval(args.starttime,args.endtime,"%Y%m%d")
TL_in = TimeList.fromfilenames(TI, INPUTDIR,args.input_glob_pattern, prefix=args.input_prefix)
TLout = TimeList.fromfilenames(TI, OUTDIR,  args.output_glob_pattern,prefix=args.output_prefix)

jpk,jpj,jpi=TheMask.shape


FILELIST_OUT=TLout.filelist[rank::nranks]
FILELIST__IN=TL_in.filelist[rank::nranks]
Timelist_OUT=TLout.Timelist[rank::nranks]



for iFrame, outfile in enumerate(FILELIST_OUT):

    inputfile = FILELIST__IN[iFrame]
    print "rank ", rank, inputfile, outfile
    VAR= DataExtractor(TheMask,inputfile,varname)
    netcdf3.write_3d_file(VAR.values, varname, outfile, TheMask, fillValue=1e+20)


