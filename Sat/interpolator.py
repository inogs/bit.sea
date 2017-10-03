import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Interpolates from 1km mesh to output mesh.
    Works in parallel
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required = True,
                                help = ''' E.g. dir with files on 1km mesh'''

                                )

    parser.add_argument(   '--outputdir', '-o',
                                type = str,
                                required = True,
                                help = ''' E.g. dir with files on 1/24 mesh'''

                                )

    parser.add_argument(   '--outmesh', '-m',
                                type = str,
                                required = True,
                                choices = ['SatOrigMesh','V4mesh','V1mesh','KD490mesh','SAT1km_mesh', 'Mesh24'],
                                help = ''' Name of the mesh of sat ORIG and used to dump checked data.'''
                                )
    parser.add_argument(   '--maskfile', '-M',
                                type = str,
                                required = True,
                                help = ''' Name of the mesh of sat ORIG and used to dump checked data.'''
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
maskOut = getattr(masks,args.outmesh)



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

x1km = Sat.OneKmMesh.lon
y1km = Sat.OneKmMesh.lat

I_START, I_END = interp2d.array_of_indices_for_slicing(x, x1km)
J_START, J_END = interp2d.array_of_indices_for_slicing(y, y1km)

INPUTDIR=addsep(args.inputdir)
OUTPUTDIR=addsep(args.outputdir)
dateformat="%Y%m%d"

reset = False

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
    Mfine = Sat.readfromfile(filename)
    Mout  = interp2d.interp_2d_by_cells_slices(Mfine, TheMask, I_START, I_END, J_START, J_END)
    Sat.dumpGenericNativefile(outfile, Mout, 'CHL', maskOut)

    print "\tfile ", counter, " of ", MySize, " done by rank ", rank