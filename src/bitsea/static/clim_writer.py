# script for Zenodo publishing
import numpy as np
from bitsea.commons.layer import Layer
import bitsea.basins.V2 as basV2
import netCDF4
from bitsea.static.climatology import get_climatology

SUBLIST = basV2.P.basin_list
PresDOWN=np.array([0,25,50,75,100,125,150,200,400,600,800,1000,1500,2000,2500])
LayerList=[ Layer(PresDOWN[k], PresDOWN[k+1])  for k in range(len(PresDOWN)-1)]

z_clim = np.array([-(l.bottom+l.top)/2  for l in LayerList])


N3n_clim, N3n_std = get_climatology('N3n', SUBLIST, LayerList,basin_expand=True)
N1p_clim, N1p_std = get_climatology('N1p', SUBLIST, LayerList,basin_expand=True)
O2o_clim, O2o_std = get_climatology('O2o', SUBLIST, LayerList,basin_expand=True)
N4n_clim, N4n_std = get_climatology('N4n', SUBLIST, LayerList,basin_expand=True)
N5s_clim, N5s_std = get_climatology('N5s', SUBLIST, LayerList,basin_expand=True)
ALK_clim, ALK_std = get_climatology('ALK', SUBLIST, LayerList,basin_expand=True)
pCO2_clim, pCO2_std = get_climatology('pCO2', SUBLIST, LayerList,basin_expand=True)
pH_clim, pH_std = get_climatology('pH', SUBLIST, LayerList,basin_expand=True)
DIC_clim, DIC_std = get_climatology('DIC', SUBLIST, LayerList,basin_expand=True)



SubDescr    =""
for i in SUBLIST : SubDescr   +=str(i.name)               + ", "

ncOUT=netCDF4.Dataset('climatologies.nc','w')
ncOUT.createDimension('subbasins',18)
ncOUT.createDimension('levels',len(z_clim))

ncvar=ncOUT.createVariable('z','f',('levels',))
ncvar[:] = -z_clim
setattr(ncvar,'units'        ,'m')
setattr(ncvar,'long_name'    ,'depth')
setattr(ncvar,'standard_name','depth')
setattr(ncvar,'positive'     ,'down')



ncvar=ncOUT.createVariable('nitrate_mean','f',('subbasins','levels'))
ncvar[:] = N3n_clim
ncvar=ncOUT.createVariable('nitrate_std','f',('subbasins','levels'))
ncvar[:] = N3n_std

ncvar=ncOUT.createVariable('phosphate_mean','f',('subbasins','levels'))
ncvar[:] = N1p_clim
ncvar=ncOUT.createVariable('phosphate_std','f',('subbasins','levels'))
ncvar[:] = N1p_std

ncvar=ncOUT.createVariable('oxygen_mean','f',('subbasins','levels'))
ncvar[:] = O2o_clim
ncvar=ncOUT.createVariable('oxygen_std','f',('subbasins','levels'))
ncvar[:] = O2o_std

ncvar=ncOUT.createVariable('silicate_mean','f',('subbasins','levels'))
ncvar[:] = N5s_clim
ncvar=ncOUT.createVariable('silicate_std','f',('subbasins','levels'))
ncvar[:] = N5s_std

ncvar=ncOUT.createVariable('alkalinity_mean','f',('subbasins','levels'))
ncvar[:] = ALK_clim
ncvar=ncOUT.createVariable('alkalinity_std','f',('subbasins','levels'))
ncvar[:] = ALK_std

ncvar=ncOUT.createVariable('partial_pressure_CO2_mean','f',('subbasins','levels'))
ncvar[:] = pCO2_clim
ncvar=ncOUT.createVariable('partial_pressure_CO2_std','f',('subbasins','levels'))
ncvar[:] = pCO2_std

ncvar=ncOUT.createVariable('pH_mean','f',('subbasins','levels'))
ncvar[:] = pH_clim
ncvar=ncOUT.createVariable('pH_std','f',('subbasins','levels'))
ncvar[:] = pH_std

ncvar=ncOUT.createVariable('DIC_mean','f',('subbasins','levels'))
ncvar[:] = DIC_clim
ncvar=ncOUT.createVariable('DIC_std','f',('subbasins','levels'))
ncvar[:] = DIC_std



setattr(ncOUT,'subbasins',SubDescr[:-2])
setattr(ncOUT,'subbasin_polygonals',"https://github.com/inogs/bit.sea/blob/master/basins/V2.py")
ncOUT.close()




