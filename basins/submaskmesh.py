from commons.mask import Mask
TheMask=Mask('/gpfs/scratch/userexternal/ateruzzi/MASKS24/meshmask.nc')
#TheMask=Mask('/gpfs/scratch/userexternal/ateruzzi/MASKS16corrected/meshmask.nc')
from basins import V0 as OGS
from commons.submask import SubMask
import numpy as np

OUTDIR = '/gpfs/scratch/userexternal/ateruzzi/MASKS24/'

already_assigned=np.zeros(TheMask.shape,dtype=np.bool)

for sub in OGS.Pred.basin_list:
    S=SubMask(sub,maskobject=TheMask)
    #already_assigned[S.mask]=True
    #ii= (S.mask) & (~already_assigned)
    S.mask[already_assigned]=False
    S.save_as_netcdf(OUTDIR + 'submask_9.nc',maskvarname=sub.name)
    already_assigned[S.mask] = True

