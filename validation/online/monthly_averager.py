import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates monthly averaged files
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help = ''' '''

                                )

    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help = ''' mask filename .'''
                                )

    parser.add_argument(   '--outdir', '-o',
                                type = str,
                                required = True,
                                help = ''' output directory'''
                                )
    return parser.parse_args()


args = argument()

from commons.Timelist import TimeInterval, TimeList
from commons.mask import Mask
from commons.time_averagers import TimeAverager3D, TimeAverager2D
import netCDF4 as NC
from commons import netcdf4
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



TL=TimeList.fromfilenames(None, INPUTDIR, "ave*N1p.nc", filtervar="N1p")


VARLIST=['P_l','O2o','N3n','P_c','N1p','ppn','pH','O3c','CO2airflux', 'pCO2','N4n','O3h','N5s','Z_c']

MONTHLY_REQS = TL.getMonthlist()
req=MONTHLY_REQS[0]
for var in VARLIST[rank::nranks]:
    indexes,weights=TL.select(req)

    #if var=='pH': var='PH'

    outfile = OUTPUTDIR + "ave." + req.string + "01-12:00:00." + var + ".nc"
    print outfile
    filelist=[]
    for k in indexes:
        t = TL.Timelist[k]
        filename = INPUTDIR + "ave." + t.strftime("%Y%m%d-%H:%M:%S") + "." + var + ".nc"
        filelist.append(filename)
    if netcdf4.dimfile(filename, var)==3:
        M3d = TimeAverager3D(filelist, weights, var, TheMask)
        netcdf4.write_3d_file(M3d, var, outfile, TheMask)
    else:
        M2d = TimeAverager2D(filelist, weights, var, TheMask)
        netcdf4.write_2d_file(M2d, var, outfile, TheMask)
