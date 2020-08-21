import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Interpolates from a fine mesh to a coarser output mesh.
    Works in parallel
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help = ''' E.g. dir with files on 1km mesh (ORIG)'''
                                )

    parser.add_argument(   '--outputdir', '-o',
                                type = str,
                                required = True,
                                help = ''' E.g. dir with files on 1/24 mesh (interpolated)'''
                                )
    parser.add_argument(   '--inmesh',
                                type = str,
                                required = True,
                                choices = ['SatOrigMesh','V4mesh','V1mesh','KD490mesh','SAT1km_mesh', 'Mesh24','Mesh4'],
                                help = ''' Name of the mesh of the input sat file'''
                                )
    parser.add_argument(   '--outmesh', '-m',
                                type = str,
                                required = True,
                                choices = ['SatOrigMesh','V4mesh','V1mesh','KD490mesh','SAT1km_mesh', 'Mesh24','Mesh4'],
                                help = ''' Name of the mesh of the output sat files'''
                                )
    parser.add_argument(   '--maskfile', '-M',
                                type = str,
                                required = True,
                                help = ''' Path of the meshmask corresponding to output sat files'''
                                )

    return parser.parse_args()

args = argument()


from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
from Sat import SatManager as Sat
from Sat import interp2d
from commons.mask import Mask
from postproc import masks
from commons.utils import addsep
import os
from commons import netcdf3
maskOut = getattr(masks,args.outmesh)
maskIn  = getattr(masks,args.inmesh)


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

TheMask = Mask(args.maskfile)

x = maskOut.lon
y = maskOut.lat

xOrig = maskIn.lon
yOrig = maskIn.lat

I_START, I_END = interp2d.array_of_indices_for_slicing(x, xOrig)
J_START, J_END = interp2d.array_of_indices_for_slicing(y, yOrig)

INPUTDIR=addsep(args.inputdir)
OUTPUTDIR=addsep(args.outputdir)
dateformat="%Y%m%d"
dateformat="%Y%m"

reset = True

Timestart="19500101"
Time__end="20500101"

TI = TimeInterval(Timestart,Time__end,"%Y%m%d")
TL = TimeList.fromfilenames(TI, INPUTDIR,"*.nc",prefix='',dateformat=dateformat)

counter = 0
MySize = len(TL.filelist[rank::nranks])

for filename in TL.filelist[rank::nranks]:
    counter += 1
    outfile = OUTPUTDIR + os.path.basename(filename)
    exit_condition = os.path.exists(outfile) and (not reset)
    if exit_condition: 
        continue
    Mfine = Sat.readfromfile(filename,'KD490')
    Mout, usedPoints  = interp2d.interp_2d_by_cells_slices(Mfine, TheMask, I_START, I_END, J_START, J_END, min_cov=0.0, ave_func=Sat.mean)
    Sat.dumpGenericNativefile(outfile, Mout, 'KD490', maskOut)
    netcdf3.write_2d_file(usedPoints, 'Points', outfile, Mout)

    print "\tfile ", counter, " of ", MySize, " done by rank ", rank

