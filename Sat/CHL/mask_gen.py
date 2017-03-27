'''
Saves the mask file for CHL 1km files
Executed once and for all.
'''

import numpy as np
from netCDF4 import Dataset
COUNT = np.load('CHL_map_occurency.npy')


filename="/gss/gss_work/DRES_OGS_BiGe/Observations/TIME_RAW_DATA/STATIC/SAT/CCI_1km/DAILY/ORIG/20061115_d-OC_CNR-L3-CHL-MedOC4AD4_SAM_1KM-MED-REP-v02.nc"

ncIN = Dataset(filename, 'r')
LON=ncIN.variables['lon'][:]
LAT=ncIN.variables['lat'][:]
ncIN.close()

tmask = COUNT > 0
jpj,jpi=COUNT.shape
ncOUT = Dataset('CHL_1km_meshmask.nc','w')
ncOUT.createDimension("lon", jpi)
ncOUT.createDimension("lat", jpj)
ncvar=ncOUT.createVariable("lon", "f", ('lon',))
ncvar[:]=LON
ncvar=ncOUT.createVariable("lat", "f", ('lat',))
ncvar[:]=LAT
ncvar=ncOUT.createVariable('tmask', 'b', ('lat','lon'))
ncvar[:]=tmask
ncOUT.close()


ncOUT = Dataset('CHL_occurrency.nc','w')
ncOUT.createDimension("lon", jpi)
ncOUT.createDimension("lat", jpj)
ncvar=ncOUT.createVariable("lon", "f", ('lon',))
ncvar[:]=LON
ncvar=ncOUT.createVariable("lat", "f", ('lat',))
ncvar[:]=LAT
ncvar=ncOUT.createVariable('nData', 'i', ('lat','lon'))
setattr(ncvar,'missing_value',0)
ncvar[:]=COUNT
ncOUT.close()




