from commons.mask import Mask
import numpy as np
import netCDF4

TheMask=Mask('/pico/home/usera07ogs/a07ogs00/OPA/V2C/etc/static-data/MED1672_cut/MASK/meshmask.nc')
CLIM_FILE="Climatology_KD490.nc"
tmask = TheMask.mask_at_level(0)

ncIN = netCDF4.Dataset(CLIM_FILE,'r')
CLIM16 = np.array( ncIN.variables['Mean'])
ncIN.close()

jpk,jpj,jpi = TheMask.shape
start_i=53

CLIM16[CLIM16==-999] = np.nan
tmask=tmask[:, start_i:]

VOIDS=np.zeros((365),np.int)


for julian in range(365):
    clim16=CLIM16[julian,:,start_i:]
    ii = (np.isnan(clim16) & tmask)
    J,I = np.nonzero(ii)
    print I
    VOIDS[julian] = ii.sum()

import sys
sys.exit()
import matplotlib.pyplot as pl
julian=0
clim16=CLIM16[julian,:,start_i:]
fig,ax = pl.subplots()
ii = (np.isnan(clim16) & tmask)
ax.imshow(tmask)
J,I =np.nonzero(ii)
for k in range(7):
    ax.plot( I[k], J[k],'w.')

ax.invert_yaxis()
fig.show()
fig.savefig("buchi_in_clim_16.png")
