from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
from commons.mask import Mask
from Sat import SatManager as Sat
from Sat.CHL import interp2d
from postproc.masks import Mesh24

import os
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

TheMask=Mask('/pico/home/userexternal/pdicerbo/WorkDir/meshmask.nc',dzvarname="e3t_0")
jpk,jpj,jpi = TheMask.shape
x = TheMask.xlevels[0,:]
y = TheMask.ylevels[:,0]

x1km = Sat.OneKmMesh.lon
y1km = Sat.OneKmMesh.lat

I_START, I_END = interp2d.array_of_indices_for_slicing(x, x1km)
J_START, J_END = interp2d.array_of_indices_for_slicing(y, y1km)

INPUTDIR="/pico/home/userexternal/pdicerbo/WorkDir/AveSat24/Checked_Weekly_Sat1km/"
OUTPUTDIR="/pico/home/userexternal/pdicerbo/WorkDir/AveSat24/Checked_Weekly_SatInterp24/"
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
    M16  = interp2d.interp_2d_by_cells_slices(Mfine, TheMask, I_START, I_END, J_START, J_END)
    Sat.dumpGenericNativefile(outfile, M16, 'CHL', Mesh24)

    print "\tfile ", counter, " of ", MySize, " done by rank ", rank