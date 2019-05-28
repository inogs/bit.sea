from commons.mask import Mask
TheMask=Mask('/gpfs/scratch/userexternal/ateruzzi/MASKS4/meshmask.nc')
from basins import OGS
from commons.submask import SubMask
import numpy as np

OUTDIR = '/gpfs/scratch/userexternal/ateruzzi/MASKS4/'

already_assigned=np.zeros(TheMask.shape,dtype=np.bool)

for sub in OGS.Pred.basin_list:
    S=SubMask(sub,maskobject=TheMask)
    #already_assigned[S.mask]=True
    #ii= (S.mask) & (~already_assigned)
    S.mask[already_assigned]=False
    S.save_as_netcdf(OUTDIR + 'submask9_4.nc',maskvarname=sub.name)
    already_assigned[S.mask] = True

