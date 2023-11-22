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
    parser.add_argument(   '--maskfile', '-m',
                                type = str,
                                required = True,
                                help = ''' Path of the meshmask corresponding to output sat files'''
                                )
    parser.add_argument(   '--varname', '-v',
                                type = str,
                                required = True,
                                choices = ['CHL','KD490','DIATO','NANO','PICO', 'DINO']
                                )
    parser.add_argument(   '--force', '-f',
                                action='store_true',
                                help = """Overwrite existing variables in files
                                """)

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
x = TheMask.xlevels[0,:]
y = TheMask.ylevels[:,0]

xOrig = maskIn.lon
yOrig = maskIn.lat

I_START, I_END = interp2d.array_of_indices_for_slicing(x, xOrig)
J_START, J_END = interp2d.array_of_indices_for_slicing(y, yOrig)

INPUTDIR=addsep(args.inputdir)
OUTPUTDIR=addsep(args.outputdir)
dateformat="%Y%m%d"

Timestart="19500101"
Time__end="20500101"

TI = TimeInterval(Timestart,Time__end,"%Y%m%d")
TL = TimeList.fromfilenames(TI, INPUTDIR,"*.nc",prefix='',dateformat=dateformat)



for filename in TL.filelist[rank::nranks]:

    outfile = OUTPUTDIR + os.path.basename(filename)
    writing_mode = Sat.writing_mode(outfile)

    condition_to_write = not Sat.exist_valid_variable(args.varname,outfile)
    if args.force: condition_to_write=True
    if not condition_to_write: continue


    Mfine = Sat.readfromfile(filename, args.varname)
    Mout, usedPoints  = interp2d.interp_2d_by_cells_slices(Mfine, TheMask, I_START, I_END, J_START, J_END, min_cov=0.0, ave_func=Sat.mean)
    Sat.dumpGenericfile(outfile, Mout, args.varname, mesh=TheMask,mode=writing_mode)
    print(outfile,flush=True)
    if not Sat.exist_valid_variable('Points', outfile):
        Sat.dumpGenericfile(outfile, usedPoints, 'Points', mesh=TheMask, mode='a')


