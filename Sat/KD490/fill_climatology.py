from commons.mask import Mask
import numpy as np
import netCDF4
import scipy.interpolate

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
CLIM16_FILLED=np.ones_like(CLIM16)*-999




for julian in range(365):
    print julian
    clim16=CLIM16[julian,:,start_i:]
    goods  = (~np.isnan(clim16) & tmask)
    Jgoods, Igoods = np.nonzero(goods)
    nP = len(Jgoods)
    points = np.zeros((nP,2),dtype=np.float32)
    points[:,0] = Jgoods
    points[:,1] = Igoods
    values = clim16[goods]
    
    tofill = (np.isnan(clim16) & tmask)
    J,I = np.nonzero(tofill)
    nP = len(J)
    xi = np.zeros((nP,2),dtype=np.float32)
    xi[:,0] = J
    xi[:,1] = I
    
    V= scipy.interpolate.griddata(points, values, xi, "nearest")
    clim16[tofill] = V
    CLIM16_FILLED[julian,:,start_i:] = clim16

CLIM16_FILLED[np.isnan(CLIM16_FILLED)] = -999.0

ncOUT = netCDF4.Dataset('Climatology_KD490_filled.nc','w')
ncOUT.createDimension('lon',jpi)
ncOUT.createDimension('lat',jpj)
ncOUT.createDimension('time',365)
ncvar=ncOUT.createVariable('Mean','f',('time','lat','lon'))
setattr(ncvar,'missing_value',-999.0)
ncvar[:]=CLIM16_FILLED
ncOUT.close()
