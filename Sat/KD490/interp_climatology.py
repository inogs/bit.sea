from commons.mask import Mask
from Sat import SatManager as Sat
import numpy as np
from Sat.KD490 import interp2d
import netCDF4

TheMask=Mask('/pico/home/usera07ogs/a07ogs00/OPA/V2C/etc/static-data/MED1672_cut/MASK/meshmask.nc')
CLIM_FILE="/gss/gss_work/DRES_OGS_BiGe/Observations/CLIMATOLOGY/SAT/KD490/Climatology_KD490.nc"
MEAN,STD = Sat.readClimatology(CLIM_FILE)

jpk,jpj,jpi = TheMask.shape
x = TheMask.xlevels[0,:]
y = TheMask.ylevels[:,0]

x1km = Sat.masks.KD490mesh.lon
y1km = Sat.masks.KD490mesh.lat

I_START, I_END = interp2d.array_of_indices_for_slicing(x, x1km)
J_START, J_END = interp2d.array_of_indices_for_slicing(y, y1km)

CLIM = np.zeros((365, jpj,jpi), dtype=[('MEAN',np.float32)])

for julian in range(365):
    print julian
    Mfine=MEAN[julian,:,:]
    M16  = interp2d.interp_2d_by_cells_slices(Mfine, TheMask, I_START, I_END, J_START, J_END)
    CLIM[julian,:,:] = M16

ncOUT = netCDF4.Dataset('Climatology_KD490.nc','w')
ncOUT.createDimension('lon',jpi)
ncOUT.createDimension('lat',jpj)
ncOUT.createDimension('time',365)
ncvar=ncOUT.createVariable('Mean','f',('time','lat','lon'))
setattr(ncvar,'missing_value',-999.0)
ncvar[:]=CLIM['MEAN']
ncOUT.close()


