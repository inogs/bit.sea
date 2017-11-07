import argparse

def argument():
    parser = argparse.ArgumentParser(description = '''
    Executes pngquant in parallel
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                required =True,
                                help = ''' Orig Maps Dir'''
                                )

    parser.add_argument(   '--outputdir', '-o',
                                type = str,
                                required = True,
                                help = 'Output compressed dir')
    parser.add_argument(   '--bindir', '-b',
                                type = str,
                                required = True,
                                help = 'HOST/ chain bin directory, where to find pngquant')

    return parser.parse_args()

args = argument()
from commons.utils import addsep
import os,glob

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

INPUTDIR =addsep(args.inputdir)
OUTPUTDIR=addsep(args.outputdir)
BINDIR   =addsep(args.bindir)

localdir = "%s%d/" %(OUTPUTDIR, rank)
os.system("mkdir -p " + localdir)

FILELIST=glob.glob(INPUTDIR + "*png")
os.chdir(localdir)
for filename in FILELIST[rank::nranks]:
    os.system("ln -fs " + filename)

os.system(BINDIR + "pngquant *png")
os.system("mv *-fs8.png ../")
os.system("rm -rf " + localdir)

    

    