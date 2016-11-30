from commons.Timelist import TimeList
from commons.time_interval import TimeInterval
from commons.mask import Mask
from Sat import SatManager as Sat
from Sat.KD490 import interp2d
import os

TheMask=Mask('/pico/home/usera07ogs/a07ogs00/OPA/V2C/etc/static-data/MED1672_cut/MASK/meshmask.nc')
jpk,jpj,jpi = TheMask.shape
x = TheMask.xlevels[0,:]
y = TheMask.ylevels[:,0]

x1km = Sat.masks.KD490mesh.lon
y1km = Sat.masks.KD490mesh.lat

I_START, I_END = interp2d.array_of_indices_for_slicing(x, x1km)
J_START, J_END = interp2d.array_of_indices_for_slicing(y, y1km)

INPUTDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/KD490/MONTHLY/ORIGMESH/"
OUTPUTDIR="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/KD490/MONTHLY/V4/"
dateformat="%Y%m"

reset = False

Timestart="19990102"
Time__end="20160101"

TI = TimeInterval(Timestart,Time__end,"%Y%m%d")
TL = TimeList.fromfilenames(TI, INPUTDIR,"*.nc",prefix='',dateformat=dateformat)

for filename in TL.filelist:
    outfile = OUTPUTDIR + os.path.basename(filename)
    exit_condition = os.path.exists(outfile) and (not reset)
    if exit_condition: 
        continue
    Mfine = Sat.readfromfile(filename,'KD490')
    M16  = interp2d.interp_2d_by_cells_slices(Mfine, TheMask, I_START, I_END, J_START, J_END)
    Sat.dump_simple_V4file(outfile, M16, 'KD490')